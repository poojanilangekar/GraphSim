
#include <mpi.h>
#include <assert.h>
#include <iostream>
#include <string>
#include <set>
#include <map>
#include <fstream>
#include <streambuf>
#include <cmath>


#include <json.hpp>

#include "graphchi_basic_includes.hpp"
using namespace graphchi;

/* 
 * query_json contains the labels and the structure of the query graph.
 */
nlohmann::json query_json;


int main(int argc, const char ** argv) {
    
    graphchi_init(argc, argv);
    
    double start_time, end_time;
    
    MPI_Init(NULL,NULL);
    
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    
    if((world_size > 1) && (world_rank == 0)) //Only one master process can exist. 
    {
        logstream(LOG_ERROR)<<"Maximum of 1 master process should be started.\n";
        exit(1);
    }

    assert(world_size == 1);
    
    start_time = MPI_Wtime();

    /*
     * Read the query file and parse the contents into a json object.
     */
    std::string queryfile = get_option_string("queryfile");
    std::ifstream qfile(queryfile);
    std::stringstream qss;
    qss << qfile.rdbuf();
    query_json = nlohmann::json::parse(qss);


    /*
     * Gather the metrics to calculate the number of machines to be used.
     * Compute the optimal number of partitions.  
     */
    int memory = get_option_int("memory",2000); //If the memory parameter is not given assume 800MB (Same as graphchi)
    if(memory <= 0)
    {
        logstream(LOG_FATAL)<<"Memory is set to"<<memory<<"mB\n";
        logstream(LOG_FATAL)<<"Memory should always be greater than 0mB\n";
        exit(1);
    }
    assert(memory>0);
    
    int hosts=0;
    std::ifstream hfile("./machines");
    std::string line;
    if(hfile.is_open()){
        while(std::getline(hfile,line))
            hosts++;
    }
    hfile.close();
    if(hosts == 0)
    {
        logstream(LOG_FATAL)<<"GraphSim requires atleast 1 machine.\n";
        exit(1);
    }
    assert(hosts>0);
    
    float alpha = 0.75;
    std::string file_type = get_option_string("filetype", std::string());
    if (file_type.empty())
        logstream(LOG_WARNING) << "The parameter \"filetype\" is not set. May cause worker process to exit if preprocessed shards are not present.\n";  
    
    size_t vertices = get_option_long("vertices",0);
    size_t edges = get_option_long("edges",0);
    
    int nparts;
    std::string worker;
    if((vertices == 0) || (edges == 0)){
        logstream(LOG_WARNING)<<"Number of edges/vertices not mentioned. Choosing default option of "<< hosts<<"workers.\n";
        nparts = hosts;
    } else {
        int qvertices = query_json["node"].size();
        int qedges = query_json["edge"].size();

        std::string vertexfile = get_option_string("vertexfile");
        std::ifstream vfile(vertexfile, std::ios::binary);
        vfile.seekg(0,vfile.end);
        size_t vfilesize = vfile.tellg();
        vfile.close();
        
        logstream(LOG_INFO)<<"Vertex file size: "<<vfilesize/1000000<<"MB\n";

        nparts = std::ceil((float)((4*edges*(32 + qvertices*4.0))/((memory*alpha*1000000) - (1152 + (vfilesize) + (vertices *( 124 + (24*qvertices)+ (4 *qedges)))) )));
        
        if((nparts < 1) || (nparts >= 10)) //If there is insufficient memory or 10 or more workers are required, use the swap worker. 
        {
            
            worker = "./bin/src/GraphSim_swap_worker";
            nparts = std::ceil(((4*edges*(32 + qvertices*4.0))+(vertices *( 120 + (24*qvertices)+ (4 *qedges)))+vfilesize)/((alpha*memory*1000000)-(1152+(vertices*4)))) + std::ceil(vertices/10000000);
            if(nparts == 2) // 2 swap workers signify insufficient memory.
                nparts = -1;
            else
                logstream(LOG_DEBUG)<<"Starting worker with swap function. Engine might be slower!\n";
        }
        else {  
            worker = "./bin/src/GraphSim_worker";
            if(nparts != 1) 
                nparts = nparts + 1; //+1 for the bloom filter.
        }
    }
    /*
     * If the number of parts calculated is less than one, the memory available is insufficient.  
     */
    if(nparts < 1)
    {
        logstream(LOG_FATAL)<<"Insufficient memory space. Try again by adding additional memory.\n";
        exit(1);
    }
    assert(nparts >= 1); 
    
    if(nparts > hosts )
    {
        logstream(LOG_FATAL)<<"Insufficient machines. Try again by adding additional machines to the cluster or increase the memory available.\n";
        exit(1);
    }
    assert(nparts <= hosts);
      

    MPI_Comm intercomm;

    MPI_Info info;

    int ierr = MPI_Info_create(&info);
    if(ierr != MPI_SUCCESS)
    {
        logstream(LOG_ERROR)<<"Error creating MPI_Info object.\n Errorcode: "<<ierr<<"\n";
        exit(1);
    }
    
    MPI_Info_set(info,"hostfile","./machines");
    
    logstream(LOG_DEBUG)<<"Spawning "<<nparts<<" worker processes.\n";

    int err = MPI_Comm_spawn(worker.c_str(),(char **) argv, nparts, info, 0, MPI_COMM_WORLD,&intercomm,NULL);

    if(err != MPI_SUCCESS)
    {
        logstream(LOG_ERROR)<<"Error creating worker processes.\n Errorcode: "<<err<<"\n";
        exit(1);
    }
    
    logstream(LOG_INFO)<<"Successfully spawned "<<nparts<<" worker processes.\n";
    
    /*
     * At the master, receive the messages from all the worker processes.
     * If the size of the result is equal to the size of the query patter, print the corresponding matches to the output file.
     * Else print 'No match' to indicate that the graph did not match the query.
     * Display the total runtime and data shipment. 
     */
    int ship_size = 0;
    std::map<vid_t, std::set <vid_t> > result;
    for (int i = 0; i < nparts; i++) {
        for (size_t j = 0; j < query_json["node"].size(); j++) {
            MPI_Status status;
            MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, intercomm, &status);
            int count;
            MPI_Get_count(&status, MPI_UINT32_T, &count);
            if ((status.MPI_SOURCE < nparts) && (status.MPI_TAG >= 0)) {
                ship_size = ship_size + count;
                vid_t *l = new vid_t[count];
                MPI_Recv(l, count, MPI_UINT32_T, status.MPI_SOURCE, status.MPI_TAG, intercomm, &status);
                logstream(LOG_INFO)<<"Received "<<count<<" matches for node "<<status.MPI_TAG<<" from worker "<<status.MPI_SOURCE<<"\n";
                if(count != 0)
                    result[status.MPI_TAG].insert(l, l + count);
                delete[] l;
            } else
                j--;
        }
    }
    std::string filename = get_option_string("file");
    std::string outputfile = get_option_string("outputfile", filename + "_output.txt");
    logstream(LOG_DEBUG)<<"Computation and Data Shipment completed! Storing the output in "<<outputfile<<"\n";
    std::ofstream output(outputfile);
    if (result.size() != query_json["node"].size())
        output << "No match";
    else {
        for (std::map <vid_t, std::set<vid_t> >::iterator it = result.begin(); it != result.end(); ++it) {
            output << it->first << "\n";
            for (std::set<vid_t>::iterator i = (it->second).begin(); i != (it->second).end(); ++i)
                output << *i << "\t";
            output << "\n";
        }
    }
    output.close();

    end_time = MPI_Wtime();
    logstream(LOG_INFO) << "Runtime: " << (end_time - start_time) << "sec\n";
    logstream(LOG_INFO) << "Data Shipped: " << sizeof (vid_t) * ship_size << "B\n";

    MPI_Finalize();
    return 0;
}

