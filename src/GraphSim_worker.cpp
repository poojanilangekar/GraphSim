//Dynamic edges. 

#define DYNAMICEDATA 1

#include <mpi.h>
#include <assert.h>
#include <iostream>
#include <string>
#include <set>
#include <map>
#include <vector>
#include <algorithm>
#include "tbb/concurrent_hash_map.h" 
#include <fstream>
#include <streambuf>
#include <thread>
#include <mutex>
#include <bloom_filter.hpp>

#include <json.hpp>

#include "graphchi_basic_includes.hpp"
#include "api/dynamicdata/chivector.hpp"

using namespace graphchi;

//Vertex Data is of int and Edge data is of type chivector<vid_t>
typedef int VertexDataType;
typedef chivector<vid_t> EdgeDataType;

/* 
 * query_json contains the labels and the structure of the query graph.
 * rvec_map is a mapping between each node and its corresponding result vectors. 
 * label_map contains the vertex labels of the input data graph.
 * bloom_schedule keeps track of all the vertices which have been explored and need to be refined.
 */

nlohmann::json query_json;
tbb::concurrent_hash_map<unsigned int, nlohmann::json> rvec_map;
tbb::concurrent_hash_map<unsigned int, nlohmann::json> label_map;
std::mutex q_mtx;
std::mutex b_mtx;
bloom_filter bloom_schedule;

void init_bloom(unsigned long int vertices) {
    bloom_parameters p;
    p.false_positive_probability = 0.75;
    p.projected_element_count = vertices;
    p.compute_optimal_parameters();
    bloom_filter temp(p);
    bloom_schedule = temp;
}

void bloom_insert(vid_t vertex_id) {
    b_mtx.lock();
    bloom_schedule.insert(vertex_id);
    b_mtx.unlock();
}

bool bloom_contains(vid_t vertex_id) {
    b_mtx.lock();
    bool contains = bloom_schedule.contains(vertex_id);
    b_mtx.unlock();
    return contains;
}

struct GraphSimulation : public GraphChiProgram<VertexDataType, EdgeDataType> {
    
/*
 * The check_equal function determines if the labels of the datanode and the querynode are equal.
 * The function iterates through all the labels of the query node and compares its value to that of the datanode.
 * It returns true only if all the labels match; returns false otherwises.
 */
     bool check_equal(nlohmann::json datanode, nlohmann::json querynode) {
         
         if(datanode.is_null())
             return false;
        for (nlohmann::json::iterator it = querynode.begin(); it != querynode.end(); ++it) {
            if ((it.key() != "id") && (it.key() != "out_degree")) {
                if (datanode.find(it.key()) != datanode.end()) {
                    if(datanode[it.key()].is_string())
                    {
                        std::string d = datanode[it.key()], q = it.value();
                        if(d.find(q)== std::string::npos)
                            return false;
                    }
                    else if (datanode[it.key()] != it.value())
                        return false;
                } else
                    return false;
            }
        }
        return true;
    }
    
    
    void update(graphchi_vertex<VertexDataType, EdgeDataType> &vertex, graphchi_context &gcontext) {
            
        /*
         * Concurrent accessor to access the rvec value corresponding to the current vertex.
         */
        tbb::concurrent_hash_map<unsigned int, nlohmann::json>::accessor ac;
        rvec_map.insert(ac,vertex.id());
        nlohmann::json rvec = ac->second;
        
        int dependencies; //The number of active dependencies of the current vertex. 
        
        /*
         * vertex_false to keep track of all the query vertices marked false for the vertex in the current iteration
         */

        std::vector<vid_t> vertex_false;

        

        /*
         * If the vertex has a null rvec, vertex is not yet explored. 
         * Compare the vertex with each of the vertices in the query graph.
         * If the current node matches the query node, add the dependencies of the query node to the rvec.
         * If the query node does not have any dependencies, set rvec[i] as true. (This implies a direct match)
         * If the query node and the current node don't match, add i to vertex_false.
         */
        
        if(rvec.is_null()){
            
            rvec = nlohmann::json::array();
            
            dependencies = 0; //Vertex is being computed for the first time and hence has zero dependencies. 
            tbb::concurrent_hash_map<unsigned int, nlohmann::json>::accessor la;
            label_map.insert(la,vertex.id());
            nlohmann::json label_data = la->second;
            la.release();
            label_map.erase(vertex.id());
            if((vertex.num_outedges()==0) && (vertex.num_inedges() == 0)) {
             vertex.set_data(0);
             ac->second = rvec;
             ac.release();
             return;
             
            }
            
            for(unsigned int i=0; i < query_json["node"].size(); i++) {
                q_mtx.lock();
                if(check_equal(label_data,query_json["node"][i])) {    
                    q_mtx.unlock();
                    unsigned int out_d = query_json["node"][i]["out_degree"];
                    if(out_d == 0){
                        rvec[i] = true;
                    }
                    else if(vertex.num_outedges() == 0)
                    {
                        vertex_false.push_back(i);
                    }
                    else
                    {   for(unsigned int j=0; j <query_json["edge"].size(); j++){
                            unsigned int source = query_json["edge"][j]["source"], target = query_json["edge"][j]["target"];
                            if(i == source )
                                rvec[i][target] = vertex.num_outedges();
                        }
                        dependencies = dependencies + out_d;
                    }

                }
                else
                {
                    q_mtx.unlock();
                    vertex_false.push_back(i);
                }
            }
            /*
             * If the vertex has dependencies, schedule the children of the current vertex (outedges).
             */
            if(dependencies != 0){
                bloom_insert(vertex.id());
                for(int i = 0; i <vertex.num_outedges();i++)
                    gcontext.scheduler->add_task(vertex.outedge(i)->vertex_id());
            }
            
            /*
             * Vertex data is set to the number of dependencies.
             * If the vertex data is greater than 0, then it is processed whenever it is scheduled in the subsequent iterations.
             * If the vertex data is 0, it is not processed in the subsequent iterations.  
             */
            vertex.set_data(dependencies);
        } 
            
        dependencies = vertex.get_data();
        
        /*
         * If the current vertex has dependencies, it has to be processed.
         * Collect the edge data of all it's outgoing edges and for each outgoing edge which is updated, update the corresponding dependency.
         * Else, clear all the outedges. 
         */
        
        if(dependencies != 0 ) {
            
            //Gather updates from all the children and update the corresponding dependencies. 
            
            nlohmann::json updates;
            
            for(int i = 0; i < vertex.num_outedges(); i++){                    
                chivector<vid_t> * e_vector = vertex.outedge(i)->get_vector();
                int evector_size = e_vector->size();
                if(bloom_contains(vertex.outedge(i)->vertex_id())) {
                    for( int j =0; j < evector_size; j++){
                        vid_t t = e_vector->get(j);
                        if(updates[t].is_null())
                            updates[t] = 1;
                        else {
                            int n = updates[t];
                            updates[t] = n +1;
                        }
                    }
                }
                e_vector->clear();
            }
            
            for(vid_t i = 0; i < updates.size(); i++ ) {
                if(updates[i].is_null())
                    continue;
                int cur_updates = updates[i];
                for(size_t j = 0; j < rvec.size(); j++){
                    if((rvec[j].is_array()))
                    {
                        if(rvec[j][i].is_number()){
                            int prev_dep = rvec[j][i];
                            if(prev_dep <= cur_updates) {
                                //rvec[j] = false; //TBR
                                rvec[j] = nlohmann::json::value_t::null;
                                vertex_false.push_back(j);
                                int out_d = query_json["node"][j]["out_degree"];
                                dependencies = dependencies - out_d;
                            }
                            else
                                rvec[j][i] = prev_dep - cur_updates;
                        }
                    }
                }
            }
            vertex.set_data(dependencies);
        } else { 
            for(int i = 0; i < vertex.num_outedges(); i++){
                chivector<vid_t> * e_vector = vertex.outedge(i)->get_vector();
                if(e_vector->size())
                    e_vector->clear();
            }
        }
   
        /*
         * If a node has been set to false in the current iteration, propagate the update through all the inedges.
         */
        if(vertex_false.size() != 0){
            for(int i=0; i < vertex.num_inedges(); i++){
                chivector<vid_t> * e_vector = vertex.inedge(i) -> get_vector();
                for(unsigned int j = 0; j < vertex_false.size(); j++)
                    e_vector->add(vertex_false[j]);
                if(bloom_contains(vertex.inedge(i)->vertex_id()) || (dependencies != 0))  //schedule parent vertices only if they are already explored or probably lie in the descendants path.
                    gcontext.scheduler->add_task(vertex.inedge(i)->vertex_id());
            }
        }
        
        // Update the result vector and release the accsessor. 
        ac->second = rvec;
        ac.release();

    }
    
    /**
     * Called before an iteration starts.
     */
    void before_iteration(int iteration, graphchi_context &gcontext) {
    }
    
    /**
     * Called after an iteration has finished.
     */
    void after_iteration(int iteration, graphchi_context &gcontext) {
      /*
       * If there were changes in the current iteration, then iterate once more.
       * If there were no changes, stop execution by setting current iteration to last iteration. 
       */
        if(gcontext.scheduler->num_tasks() != 0)
            gcontext.set_last_iteration(iteration+1);
        else {
            gcontext.set_last_iteration(iteration);
            label_map.clear();
        }
    }
    
    /**
     * Called before an execution interval is started.
     */
    void before_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &gcontext) {        
    }
    
    /**
     * Called after an execution interval has finished.
     */
    void after_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &gcontext) {        
    }
    
};

/*
 * Function to fill the outdegree of the query nodes. 
 */
void fill_degree(nlohmann::json &query)
{
    for( size_t i = 0; i < query["node"].size(); i++) {
        query["node"][i]["out_degree"] = 0;
    }
    for( size_t i = 0; i < query["edge"].size(); i++)
    {
        size_t snode = query["edge"][i]["source"];
        query["node"][snode]["out_degree"] = int(query["node"][snode]["out_degree"]) + 1;
        
    }
}    

/*
 * Function to parse the vertexfile.
 * The label data of the input graph is stored in a json object. 
 */
void fill_vertex(std::string vfilename)
{
    std::ifstream vf(vfilename);
    std::stringstream vss;
    vss << vf.rdbuf();
    nlohmann::json vertex_json = nlohmann::json::parse(vss);
    tbb::concurrent_hash_map<unsigned int, nlohmann::json>::accessor ac;
    for(unsigned int i = 0; i <vertex_json.size(); i++) {
        label_map.insert(ac, i);
        ac->second = vertex_json[i];
        ac.release();
    }
}

int main(int argc, const char ** argv) {
    
    //Initialize worker process
    graphchi_init(argc, argv);
    int provided;    
    MPI_Init_thread(NULL,NULL, MPI_THREAD_FUNNELED, &provided);
    if(provided < MPI_THREAD_FUNNELED)
    {
        std::cout<<"Error: The MPI Library does not have support multithreaded processes.\n";
        exit(1);
    }
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int pname_len;
    
    MPI_Get_processor_name(processor_name,&pname_len);
    logstream(LOG_INFO) << "Launched worker with id "<<world_rank<<" at node "<<processor_name<<"\n";
    
    
    /*
     * Redirect stdout and stderr to the log file identified by the rank of the worker process. 
     */
    FILE* f = freopen(std::string("LOG_"+std::to_string(world_rank)).c_str(),"w",stdout);
    assert(f != NULL);
    dup2(fileno(stdout), fileno(stderr));
    
    MPI_Comm master;
    MPI_Comm_get_parent(&master); //Communicator of Master, to ship results back to the master.
    
    std::string vertexfile = get_option_string("vertexfile"); //Vertex file
    std::thread thread_fill_vertex(fill_vertex, vertexfile); 
    
    std::string queryfile = get_option_string("queryfile"); //Query file
    std::ifstream qfile(queryfile);
    std::stringstream qss;
    qss << qfile.rdbuf();
    query_json = nlohmann::json::parse(qss);
    
    std::string filename = get_option_string("file");  // Base filename
    metrics m("graph_simulation");


    bool scheduler       = true; // Always enable scheduling.

    //Shard creation.
    std::string file_type = get_option_string("filetype", std::string());
    int nshards          = convert_if_notexists<vid_t>(filename, std::to_string(world_size),file_type);
    
    /* 
     * Parse the queryfile.
     * The queryfile contains the query pattern. It includes the label data and the structure of the query.
     */
    fill_degree(query_json);
    int niters = query_json["node"].size()*query_json["edge"].size(); 
    
    
    //Initailize the graph simulation program. 
    GraphSimulation program;
    graphchi_engine<VertexDataType, EdgeDataType> engine(filename, nshards, scheduler, m); 
    

    /*
     * Optimize the graphchi engine and run the graph_simulation program.
     */
    int n_threads = std::thread::hardware_concurrency();

    engine.set_exec_threads(n_threads);
    engine.set_load_threads(n_threads);
    int memory = std::ceil((float)(get_option_int("memory",1024) * get_option_float("alpha",0.75)));
    engine.set_membudget_mb(memory);
    init_bloom(engine.num_vertices()/nshards);
    thread_fill_vertex.join();
    engine.set_reset_vertexdata(true);

    if(world_size > 1)
        engine.run(program,niters,world_rank);
    else
        engine.run(program,niters); //If world_size is 1, the engine is run using the inmemorymode. Need not mention starting interval. 


    /* After completion of the graphchi program, ship the results to the master. 
     * Construct a map result which contains a mapping between each node in the query pattern and the possible matches from the input graph.
     * If the size of the query pattern is equal to the size of the result match, then query is found. Ship the result arrays with the corresponding tags. 
     * Else the data graph does not match the query. Ship an empty message to the master.
     */
    
    
    std::map<vid_t, std::set<vid_t> > result;
    for (tbb::concurrent_hash_map<vid_t, nlohmann::json>::const_iterator it = rvec_map.begin(); it != rvec_map.end() ; ++it) {
        nlohmann::json rvec = it->second;
        for(unsigned int i = 0; i < rvec.size(); i++){
            if((rvec[i].is_array()) || (rvec[i].is_boolean())){
                result[i].insert(it->first);
            }
        }
   }
    if(result.size() != query_json["node"].size())
    {
        vid_t *p_array = NULL;
        for(size_t i = 0; i < query_json["node"].size(); i++)
            MPI_Send(p_array,0,MPI_UINT32_T,0,i,master);
    }
    else {

        for(size_t i = 0; i < query_json["node"].size();i++){
            vid_t *p_array = new vid_t[result[i].size()];
            std::copy(result[i].begin(),result[i].end(),p_array);
            MPI_Send(p_array,result[i].size(),MPI_UINT32_T,0,i,master);
            delete[] p_array;
        }
    }
    metrics_report(m);
    engine.reinitialize_edge_data(0);
    
    MPI_Finalize();
    return 0;
}
