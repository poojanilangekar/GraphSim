#GraphSim

Large-scale graph processing is becoming central to our modern life. For instance, graph pattern matching  (GPM)  can  be  utilized  to  search  and  analyze  social  graphs,  biological  data  and  road networks,  to  mention  a  few.  Conceptually, a GPM  algorithm is typically  defined  in  terms  of subgraph isomorphism, whereby it seeks to find subgraphs in an input data graph, *G*, which are similar to a given query graph, *Q*. Although subgraph isomorphism forms a uniquely important class  of  graph queries, it  is NP-complete  and  very  restrictive. Consequently, GPM has been relaxed and defined in terms of **Graph Simulation**. As opposed to subgraph isomorphism, graph simulation can run in quadratic time, return more intuitive matches, and scale well with modern big graphs(i.e., graphs with billions of verticesand edges)

GraphSim is an adaptive graph simulation system, which can run with different numbers of machines for different data and query graphs (i.e., workloads). GraphSim accomplishes its goals via employing new data model that **entirely avoids intermediate data shipment** between slave machines. It deploys a new computation model that **increases parallelism and improves memory utilization**. It adopts a new mathematical model that allows **predicting different numbers of machines** for workloads with varying complexities.

GraphSim is implemented in C++, and available as open-source under the MIT License.

Downloading GraphSim
-------------

You can download GraphSim directly from the Github Repository. Github also offers a zip download of the repository if you do not have git.

The git command line for cloning the repository is:

git clone https://github.com/poojanilangekar/GraphSim.git
cd graphlab



Building
------------------
The current version of GraphSim was tested on Ubuntu Linux 64-bit 14.04. It requires a 64-bit operating system. 
 

Dependencies
------------------

GraphSim has the following dependencies.

1. [g++ (>= 4.8)](https://gcc.gnu.org/gcc-4.8/)

2. [MPICH (>= 3.1)](https://www.mpich.org/downloads/)

3. [zlib](https://launchpad.net/ubuntu/+source/zlib)

4. [TBB](https://www.threadingbuildingblocks.org/) 


Usage 
----------------
The //Wiki-entry provides a through guide to install and run GraphSim on a cluster.


Todo
---------

The following extensions will be made to GraphSim  

	1. Solve other notions of Graph Simulation.	
	2. Enable Data Base search for label checking.  
	3. Develop a job scheduler for GraphSim.

Contributing
-------------------
1. Fork it ( https://github.com/[my-github-username]/GraphSim/fork )
2. Create your feature branch (`git checkout -b my-new-feature`)
3. Commit your changes (`git commit -am 'Add some feature'`)
4. Push to the branch (`git push origin my-new-feature`)
5. Create a new Pull Request
