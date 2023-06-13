//
// Created by milinda on 9/3/16.
//

/**
* @author Milinda Shayamal Fernando
* School of Computing , University of Utah
*
* @breif Sequential program to evaluate the performance of the local sort and SFC::seq::bucketing function.
* @details This was written to debug the weak sacalability issues with SFC::treeSort and bucketing function.
*/
#include "bucketBench.h"

int main(int  argc, char** argv)
{


    MPI_Init(&argc, &argv);

    int rank, npes;
    MPI_Comm GLOBAL_COMM=MPI_COMM_WORLD;
    MPI_Comm_rank(GLOBAL_COMM, &rank);
    MPI_Comm_size(GLOBAL_COMM, &npes);

    unsigned int options=0;


    if(argc<3)
    {   if(!rank)
            std::cout<<"Usage :"<<argv[0]<<" numPts "<<" dim "<<" maxDepth "<<" distribution (0- Normal 1- Uniform 2- LogarithmicNormal)"<<std::endl;
    }

    DendroIntL numPts=atoll(argv[1]);
    unsigned int dim=atoi(argv[2]);
    unsigned int maxDepth=atoi(argv[3]);
    unsigned int distribution=0; // default set to the gaussian distribution
    if(argc>4)
        distribution=atoi(argv[4]);

    _InitializeHcurve(dim);
    if(!rank) std::cout << "Initialized Hcurves for dimention "<<dim << std::endl;

    DendroIntL totPts=numPts*dim;
    double * pts=new double[totPts];
    std::vector<ot::TreeNode> tmpNodes;

    /*std::vector<double> pts;

    char FileName[256];
    sprintf(FileName, "%s%d_%d.pts", "ip", rank, npes);

    genGauss(0.1,numPts,dim,"ip",GLOBAL_COMM);*/

    if(distribution==0)
        genGauss(0.1,numPts,dim,pts);
    else if(distribution==1)
        genUniformRealDis(0.1,numPts,dim,pts);
    else if(distribution==2)
        genLogarithmicGauss(0.5,numPts,dim,pts);
    else
        genGauss(0.5,numPts,dim,pts);  // Default case.


    pts2Octants(tmpNodes,pts,totPts,dim,maxDepth); // generate octants from the points.
    if (!rank) std::cout << "finished generating "<<dim<<" dimensional points" << std::endl;

    std::vector<ot::TreeNode> tmpNodes_cpy;
    tmpNodes_cpy=tmpNodes;

    std::vector<ot::TreeNode> tmpSorted;
    std::vector<ot::TreeNode> tmpConstruct;
    std::vector<ot::TreeNode> tmpBalance;

    ot::TreeNode root(0,0,0,0,dim,maxDepth);

    auto t1=std::chrono::high_resolution_clock::now();
    SFC::seqSort::SFC_treeSort((&(*(tmpNodes.begin()))),tmpNodes.size(),tmpSorted,tmpConstruct,tmpBalance,maxDepth,maxDepth,root,0,1,1);
    double local_sort =std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t1).count();//MPI_Wtime()-t1;

    double stat_property[3];

    par::Mpi_Reduce(&local_sort,&stat_property[0],1,MPI_MIN,0,GLOBAL_COMM);
    par::Mpi_Reduce(&local_sort,&stat_property[1],1,MPI_SUM,0,GLOBAL_COMM);
    par::Mpi_Reduce(&local_sort,&stat_property[2],1,MPI_MAX,0,GLOBAL_COMM);



    if(!rank)
    {
        stat_property[1]=stat_property[1]/npes;
        std::cout<<"npes\tls_min\tls_mean\tls_max"<<std::endl;
        std::cout<<npes<<"\t"<<stat_property[0]<<"\t"<<stat_property[1]<<"\t"<<stat_property[2]<<std::endl;
    }


    DendroIntL splitterTemp[(NUM_CHILDREN+1)];
    tmpNodes=tmpNodes_cpy;

    DendroIntL begin=0;
    DendroIntL end=tmpNodes.size();

    t1=std::chrono::high_resolution_clock::now();
    SFC::seqSort::SFC_bucketing(&(*(tmpNodes.begin())),0,maxDepth,0,begin,end,splitterTemp);
    double bucket_time=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t1).count();//MPI_Wtime()-t1;

    par::Mpi_Reduce(&bucket_time,&stat_property[0],1,MPI_MIN,0,GLOBAL_COMM);
    par::Mpi_Reduce(&bucket_time,&stat_property[1],1,MPI_SUM,0,GLOBAL_COMM);
    par::Mpi_Reduce(&bucket_time,&stat_property[2],1,MPI_MAX,0,GLOBAL_COMM);



    if(!rank)
    {
        stat_property[1]=stat_property[1]/npes;
        std::cout<<"npes\tbkt_min\tbkt_mean\tbkt_max"<<std::endl;
        std::cout<<npes<<"\t"<<stat_property[0]<<"\t"<<stat_property[1]<<"\t"<<stat_property[2]<<std::endl;
    }

    unsigned int sf_k=16;
    tmpNodes=tmpNodes_cpy;
    t1=std::chrono::high_resolution_clock::now();
    SFC::parSort::SFC_treeSort(tmpNodes,tmpSorted,tmpConstruct,tmpBalance,0.1,maxDepth,root,0,1,1,sf_k,GLOBAL_COMM);
    double ts_time =std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t1).count();//MPI_Wtime()-t1;

    par::Mpi_Reduce(&ts_time,&stat_property[0],1,MPI_MIN,0,GLOBAL_COMM);
    par::Mpi_Reduce(&ts_time,&stat_property[1],1,MPI_SUM,0,GLOBAL_COMM);
    par::Mpi_Reduce(&ts_time,&stat_property[2],1,MPI_MAX,0,GLOBAL_COMM);



    if(!rank)
    {
        stat_property[1]=stat_property[1]/npes;
        std::cout<<"npes\tts_min\tts_mean\tts_max"<<std::endl;
        std::cout<<npes<<"\t"<<stat_property[0]<<"\t"<<stat_property[1]<<"\t"<<stat_property[2]<<std::endl;
    }


    MPI_Finalize();

    return 0;
}