//
// Created by milinda on 8/8/16.
//
#include "octUtilsBenchmark.h"





void sfcTreeSortTest(DendroIntL numPts, unsigned int dim, unsigned int maxDepth, double tolerance, int distribution,
                     unsigned int k, unsigned int options,unsigned int sf_k, char *prefix, MPI_Comm pcomm)
{


    int rank,npes;

    MPI_Comm_rank(pcomm,&rank);
    MPI_Comm_size(pcomm,&npes);

    char profFileName[256];
    sprintf(profFileName,"%s%d_%d.stat",prefix,rank,npes);
    std::vector<ot::TreeNode> tmpNodes;
    DendroIntL totPts=numPts*dim;
    double * pts=new double[totPts];

    if(distribution==0)
        genGauss(0.15,numPts,dim,pts);
    else if(distribution==1)
        genUniformRealDis(0.15,numPts,dim,pts);
    else if(distribution==2)
        genLogarithmicGauss(0.15,numPts,dim,pts);
    else
        genGauss(0.15,numPts,dim,pts);  // Default case.

    pts2Octants(tmpNodes,pts,totPts,dim,maxDepth); // generate octants from the points.

    std::vector<ot::TreeNode> pNodesSorted;
    std::vector<ot::TreeNode> pNodesConstructed;
    std::vector<ot::TreeNode> pNodesBalanced;

    pNodesSorted.clear();
    pNodesConstructed.clear();
    pNodesBalanced.clear();
    ot::TreeNode root(0,0,0,0,dim,maxDepth);

    SFC::parSort::SFC_treeSort(tmpNodes,pNodesSorted,pNodesConstructed,pNodesBalanced,tolerance,maxDepth,root,ROOT_ROT_ID,k,options,sf_k,pcomm);
    treeNodesTovtk(pNodesConstructed,rank,"oct_0.0001");
    tmpNodes.clear();


}




void weakScalingDriver(DendroIntL numPts,unsigned int dim,unsigned int maxDepth,double tolerance,int distribution,unsigned int k,unsigned int options,unsigned int sf_k,char * prefix,MPI_Comm comm)
{

    int rank,npes;

    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&npes);


#ifdef PROFILE_TREE_SORT
    unsigned int proc_group_index=binOp::fastLog2(npes);
    RuntimeStat weakScalingResults[proc_group_index]; // contains all the weak scaling results.
#endif


    double t1=0;
    double warmup_time=0;
    if(!rank) std::cout<<"======================================================================"<<std::endl;
    if(!rank) std::cout<<"     Warm Up run  (par) begin                               "<<std::endl;
    if(!rank) std::cout<<"======================================================================"<<std::endl;

    //t1=MPI_Wtime();
    sfcTreeSortTest(numPts,dim,maxDepth,tolerance,distribution,k,options,sf_k,prefix,comm);
#ifdef PROFILE_TREE_SORT
    stats.clear();
    stats_previous.clear();
#endif

   /* warmup_time=MPI_Wtime()-t1;
    std::cout<<"rank: "<<rank<<"warmup time: "<<warmup_time<<std::endl;*/


    if(!rank) std::cout<<"======================================================================"<<std::endl;
    if(!rank) std::cout<<"     Warm Up run  (par) end                               "<<std::endl;
    if(!rank) std::cout<<"======================================================================"<<std::endl;

    MPI_Barrier(comm);


#ifdef RUN_WEAK_SCALING

    int proc_group = 0;

    int min_np = 2;
    for (int i = npes; rank < i && i >= min_np; i = i >> 1) proc_group++;
    MPI_Comm comm_ws;

    MPI_Comm_split(comm, proc_group, rank, &comm_ws);

    MPI_Comm_rank(comm_ws, &rank);
    MPI_Comm_size(comm_ws, &npes);
#ifdef PROFILE_TREE_SORT
    proc_group_index=binOp::fastLog2(npes)-1;

    if(!rank) {

        weakScalingResults[proc_group_index].npes = npes;
        weakScalingResults[proc_group_index].numPts = numPts;
        weakScalingResults[proc_group_index].options = options;
        weakScalingResults[proc_group_index].tolerance = tolerance;
    }
#endif
    MPI_Barrier(comm_ws);
    /*if(!rank) std::cout<<"======================================================================"<<std::endl;
    if(!rank) std::cout<<"     TreeSort Test  (par) begin                               " <<npes<<std::endl;
    if(!rank) std::cout<<"======================================================================"<<std::endl;*/

    sfcTreeSortTest(numPts,dim,maxDepth,tolerance,distribution,k,options,sf_k,prefix,comm_ws);
#ifdef PROFILE_TREE_SORT
    if(!rank) {
        weakScalingResults[proc_group_index].stats_prev = stats_previous;
        weakScalingResults[proc_group_index].stats = stats;
    }
#endif
    /*if(!rank) std::cout<<"======================================================================"<<std::endl;
    if(!rank) std::cout<<"     TreeSort Test  (par) end                               " <<npes<<std::endl;
    if(!rank) std::cout<<"======================================================================"<<std::endl;*/

#endif

    // full run
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    if(!rank) std::cout<<"======================================================================"<<std::endl;
    if(!rank) std::cout<<"     Weak Scale Run Complete. Beginning Full Run.   "<<std::endl;
    if(!rank) std::cout<<"======================================================================"<<std::endl;

    MPI_Barrier(comm);

#ifdef PROFILE_TREE_SORT
    proc_group_index=binOp::fastLog2(npes)-1;

    if(!rank) {
        weakScalingResults[proc_group_index].npes = npes;
        weakScalingResults[proc_group_index].numPts = numPts;
        weakScalingResults[proc_group_index].options = options;
        weakScalingResults[proc_group_index].tolerance = tolerance;
    }
#endif
    /*if(!rank) std::cout<<"======================================================================"<<std::endl;
    if(!rank) std::cout<<"     TreeSort Test  (par) begin                               " <<npes<<std::endl;
    if(!rank) std::cout<<"======================================================================"<<std::endl;*/

    sfcTreeSortTest(numPts,dim,maxDepth,tolerance,distribution,k,options,sf_k,prefix,comm);

#ifdef PROFILE_TREE_SORT
    if(!rank) {
        weakScalingResults[proc_group_index].stats_prev = stats_previous;
        weakScalingResults[proc_group_index].stats = stats;
    }
#endif
    /*if(!rank) std::cout<<"======================================================================"<<std::endl;
    if(!rank) std::cout<<"     TreeSort Test  (par) end                               " <<npes<<std::endl;
    if(!rank) std::cout<<"======================================================================"<<std::endl;*/

    if(!rank) std::cout<<"======================================================================"<<std::endl;
    if(!rank) std::cout<<"     Benchmark Execution Complete.   "<<std::endl;
    if(!rank) std::cout<<"======================================================================"<<std::endl;



#ifdef PROFILE_TREE_SORT
    proc_group_index=binOp::fastLog2(npes);
    char fileName[256];
    std::ofstream statFile;
     for(unsigned int i=0;i<proc_group_index;i++)
     {

            if((options & TS_REMOVE_DUPLICATES) | options==0){
                int SF_Stages=0;
                int SF_stage_begin=0;
                if(weakScalingResults[i].npes) {

                    sprintf(fileName, "%s_%d.stat", prefix, weakScalingResults[i].npes);
                    statFile.open(fileName);
                    SF_Stages=std::ceil((binOp::fastLog2(npes)/(double)binOp::fastLog2(sf_k))) - 1;

                    if (weakScalingResults[i].npes == 2) {
                        statFile<< "npes\tnumPts\toption\ttol\t1_sf_min\t1_sf_mean\t1_sf_max\t1_split_min\t1_split_mean\t1_split_max\t1_a2a1_min\t1_a2a1_mean\t1_a2a1_max\t1_a2a2_min\t1_a2a2_mean\t1_a2a2_max\t1_ls_min\t1_ls_mean\t1_ls_max\t1_sq::rd_min\t1_sq::rd_mean\t1_sq::rd_max\t1_auxBal_min\t1_auxBal_mean\t1_auxBal_max\t1_par::rd_min\t1_par::rd_mean\t1_par::rd_max\t1_totRD_min\t1_totRD_mean\t1_totRD_max";

                        for(int i=0;i<SF_Stages;i++)
                            statFile<<"\tsf"<<(i+1)<<"_full_min\t"<<"sf"<<(i+1)<<"_full_mean\t"<<"sf"<<(i+1)<<"_full_max\t"<<"sf"<<(i+1)<<"_all2all_min\t"<<"sf"<<(i+1)<<"_all2all_mean\t""sf"<<(i+1)<<"_all2all_max\t"<<"sf"<<(i+1)<<"_split_min\t"<<"sf"<<(i+1)<<"_split_mean\t"<<"sf"<<(i+1)<<"_split_max";

                        statFile<<std::endl;



                    }
                    statFile << weakScalingResults[i].npes << "\t" << weakScalingResults[i].numPts << "\t" <<
                    weakScalingResults[i].options << "\t" << weakScalingResults[i].tolerance;
                    for (unsigned int j = 0; j < weakScalingResults[i].stats.size(); j++)
                        statFile << "\t" << weakScalingResults[i].stats[j];
                    SF_stage_begin=std::ceil((binOp::fastLog2(weakScalingResults[i].npes)/(double)binOp::fastLog2(sf_k))) - 1;
                    if(SF_stage_begin<0)
                        SF_stage_begin=0;
                    for(int k=SF_stage_begin;k<SF_Stages;k++)
                        statFile<<"\t0\t0\t0\t0\t0\t0\t0\t0\t0";
                    statFile << std::endl;
                    statFile.close();
                }


            }else
            {
                if(weakScalingResults[i].npes) {

                    sprintf(fileName,"%s_%d.stat",prefix,weakScalingResults[i].npes);
                    statFile.open (fileName);
                    if(weakScalingResults[i].npes==2)
                        statFile<<"npes\tnumPts\toption\ttol\t1_sf_min\t1_sf_mean\t1_sf_max\t1_split_min\t1_split_mean\t1_split_max\t1_a2a1_min\t1_a2a1_mean\t1_a2a1_max\t1_a2a2_min\t1_a2a2_mean\t1_a2a2_max\t1_ls_min\t1_ls_mean\t1_ls_max\t1_sq::rd_min\t1_sq::rd_mean\t1_sq::rd_max\t1_auxBal_min\t1_auxBal_mean\t1_auxBal_max\t1_par::rd_min\t1_par::rd_mean\t1_par::rd_max\t1_totRD_min\t1_totRD_mean\t1_totRD_max\t2_sf_min\t2_sf_mean\t2_sf_max\t2_split_min\t2_split_mean\t2_split_max\t2_a2a1_min\t2_a2a1_mean\t2_a2a1_max\t2_a2a2_min\t2_a2a2_mean\t2_a2a2_max\t2_ls_min\t2_ls_mean\t2_ls_max\t2_sq::rd_min\t2_sq::rd_mean\t2_sq::rd_max\t2_auxBal_min\t2_auxBal_mean\t2_auxBal_max\t2_par::rd_min\t2_par::rd_mean\t2_par::rd_max\t2_totRD_min\t2_totRD_mean\t2_totRD_max\t2_tot_min\t2_tot_mean\t2_tot_max"<<std::endl;

                    statFile<<weakScalingResults[i].npes<<"\t"<<weakScalingResults[i].numPts<<"\t"<<weakScalingResults[i].options<<"\t"<<weakScalingResults[i].tolerance;
                    for(unsigned int j=0;j<weakScalingResults[i].stats_prev.size();j++)
                        statFile<<"\t"<<weakScalingResults[i].stats_prev[j];
                    for(unsigned int j=0;j<weakScalingResults[i].stats.size();j++)
                        statFile<<"\t"<<weakScalingResults[i].stats[j];
                    statFile<<std::endl;
                    statFile.close();

                }

            }

     }
#endif


}


int main(int argc, char** argv) {

    MPI_Init(&argc, &argv);
    int rank, npes;
    MPI_Comm GLOBAL_COMM=MPI_COMM_WORLD;
    MPI_Comm_rank(GLOBAL_COMM, &rank);
    MPI_Comm_size(GLOBAL_COMM, &npes);

    unsigned int options=0;


    if(argc<3)
    {   if(!rank)
            std::cout<<"Usage :"<<argv[0]<<" numPts "<<" dim "<<" maxDepth "<<" tol  "<<" distribution (0- Normal 1- Uniform 2- LogarithmicNormal) options( 1-Remove duplicates 2- constructOctree 6-balancedOctree) "<<"statPrefix SplitterFixK"<<std::endl;
    }

    DendroIntL numPts=atoll(argv[1]);
    unsigned int dim=atoi(argv[2]);
    unsigned int maxDepth=atoi(argv[3]);
    double tol=0.001;
    unsigned int distribution=0;
    unsigned int sf_k=2;
    if(argc>4)
        tol=atof(argv[4]);

    if(argc>5)
        distribution=atoi(argv[5]);

    if(argc>6)
        options=atoi(argv[6]);

    if(argc>8)
        sf_k=atoi(argv[8]);

    if(!rank) std::cout<<"sf parameter: "<<sf_k<<std::endl;

    _InitializeHcurve(dim);
    if(!rank) std::cout << "Initialized H-Curves for dimension "<<dim << std::endl;
    weakScalingDriver(numPts,dim,maxDepth,tol,distribution,1,options,sf_k,argv[7],MPI_COMM_WORLD);




    MPI_Finalize();



}