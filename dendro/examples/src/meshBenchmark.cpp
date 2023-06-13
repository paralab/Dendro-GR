//
// Created by milinda on 9/8/16.
//

#include "meshBenchmark.h"



void meshBenchMark(char * ptsFile,bool genPts,unsigned int numPts,unsigned int dim, unsigned int maxDepth,unsigned int distribution,double tol,unsigned int sf_k,unsigned int options,char * prefix, MPI_Comm comm)
{

    int rank,npes;
    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&npes);

    std::vector<ot::TreeNode> pNodesSorted;
    std::vector<ot::TreeNode> pNodesConstructed;
    std::vector<ot::TreeNode> pNodesBalanced;

    std::vector<ot::TreeNode> tmpNodes;
    m_uiMaxDepth=maxDepth; // this is needed because we increase this when we embedded the tree.

    double t_rd=0;
    double t_cons=0;
    double t_bal=0;
    double t_mesh=0;


    double t_rd_g[3];
    double t_cons_g[3];
    double t_bal_g[3]; //0 -min 1-mean 2-Max
    double t_mesh_g[3]; //0 -min 1-mean 2-Max

    DendroIntL numUniqueOct;
    DendroIntL numConsOct;
    DendroIntL numBalOct;
    DendroIntL localSz;
    DendroIntL globalSz;
    DendroIntL cg_sz;




    unsigned int stencilSz=1;
    unsigned int eleOrder=4;

    DendroIntL totPts=numPts*dim;
    double * pts=new double[totPts];
    std::vector<double> ptsFromFile;
    char pFile[256];
    if(genPts) {

        if (distribution == 0)
            genGauss(0.15,numPts,dim,ptsFile,comm);
        else if (distribution == 1)
            genUniformRealDis(0.1, numPts, dim, pts);
        else if (distribution == 2)
            genLogarithmicGauss(0.5, numPts, dim, pts);
        else
            genGauss(0.5, numPts, dim, pts);  // Default case.
    }

    sprintf(pFile, "%s%d_%d.pts", ptsFile, rank, npes);
    //std::cout<<"Attempt to Read "<<ptsFileName<<std::endl;

    //Read pts from files
    if (!rank) {
        std::cout << RED " Reading  " << ptsFile << NRM << std::endl; // Point size
    }
    IO::readPtsFromFile(pFile, ptsFromFile);

    if (!rank) {
        std::cout << GRN " Finished reading  " << ptsFile << NRM << std::endl; // Point size
    }
    pts2Octants(tmpNodes,&(*(ptsFromFile.begin())),totPts,dim,maxDepth); // generate octants from the points.
    ptsFromFile.clear();
    delete[] pts;
    pts=NULL;

    if(!rank) {
        std::cout << YLW << "grain Size: " << numPts << std::endl;
        std::cout<<"npes: "<<npes<<std::endl;
        std::cout<<"dim: "<<dim<<std::endl;
        std::cout<<"maxDepth: "<<maxDepth<<std::endl;
        std::cout<<"tolerance: "<<tol<<std::endl;
        std::cout<<"splitter_fix_k: "<<sf_k<<std::endl;
        std::cout<<"stencil_size: "<<stencilSz<<std::endl;
        std::cout<<"distribution: "<<distribution<<NRM<<std::endl;

    }

    pNodesSorted.clear();
    pNodesConstructed.clear();
    pNodesBalanced.clear();
    ot::TreeNode root(0,0,0,0,dim,maxDepth);

    if(npes>=2) {


        if (!rank) std::cout << RED << "Remove duplicates begin" << NRM << std::endl;
        auto t1=MPI_Wtime();
        SFC::parSort::SFC_treeSort(tmpNodes, pNodesSorted, pNodesConstructed, pNodesBalanced, tol, maxDepth, root,ROOT_ROT_ID, 1, TS_REMOVE_DUPLICATES, sf_k, comm);
        auto t2=MPI_Wtime();
        t_rd=t2-t1;
        if (!rank) std::cout << RED << "Remove duplicates end" << NRM << std::endl;
        tmpNodes.clear();

        par::Mpi_Reduce(&t_rd, t_rd_g, 1, MPI_MIN, 0, comm);
        par::Mpi_Reduce(&t_rd, t_rd_g + 1, 1, MPI_SUM, 0, comm);
        par::Mpi_Reduce(&t_rd, t_rd_g + 2, 1, MPI_MAX, 0, comm);
        t_rd_g[1] = t_rd_g[1] / npes;

        localSz=pNodesSorted.size();
        par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);
        numUniqueOct=globalSz;
        if(!rank)
        {
            std::cout << YLW<< "Number of unique octants: "<<numUniqueOct<<NRM<<std::endl;
            std::cout << YLW<< "Remove duplicates time(max): "<<t_rd_g[2]<<NRM<<std::endl;
        }

        MPI_Barrier(comm);

        if (!rank) std::cout << RED << "Tree construction begin" << NRM << std::endl;
        t1=MPI_Wtime();
        SFC::parSort::SFC_treeSort(pNodesSorted, pNodesSorted, pNodesConstructed, pNodesBalanced, tol, maxDepth, root,ROOT_ROT_ID, 1, TS_CONSTRUCT_OCTREE, sf_k, comm);
        t2=MPI_Wtime();
        t_cons=t2-t1;
        if (!rank) std::cout << RED << "Tree construction end" << NRM << std::endl;
        pNodesSorted.clear();

        par::Mpi_Reduce(&t_cons, t_cons_g, 1, MPI_MIN, 0, comm);
        par::Mpi_Reduce(&t_cons, t_cons_g + 1, 1, MPI_SUM, 0, comm);
        par::Mpi_Reduce(&t_cons, t_cons_g + 2, 1, MPI_MAX, 0, comm);
        t_cons_g[1] = t_cons_g[1] / npes;


        localSz=pNodesConstructed.size();
        par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);
        numConsOct=globalSz;
        if(!rank)
        {
            std::cout << YLW<< "Number of unbalanced octants: "<<numConsOct<<NRM<<std::endl;
            std::cout << YLW<< "Tree construction time(max): "<<t_cons_g[2]<<NRM<<std::endl;
        }

        MPI_Barrier(comm);

        if (!rank) std::cout << RED << "2:1 balance begin" << NRM << std::endl;
        t1 = MPI_Wtime();//std::chrono::high_resolution_clock::now();
        SFC::parSort::SFC_treeSort(pNodesConstructed, pNodesSorted, pNodesConstructed, pNodesBalanced, tol, maxDepth, root,ROOT_ROT_ID, 1, TS_BALANCE_OCTREE, sf_k, comm);
        t2 = MPI_Wtime();//std::chrono::high_resolution_clock::now();
        t_bal = t2-t1;//std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        if (!rank) std::cout << RED << "2:1 balance end" << NRM << std::endl;
        pNodesConstructed.clear();

        par::Mpi_Reduce(&t_bal, t_bal_g, 1, MPI_MIN, 0, comm);
        par::Mpi_Reduce(&t_bal, t_bal_g + 1, 1, MPI_SUM, 0, comm);
        par::Mpi_Reduce(&t_bal, t_bal_g + 2, 1, MPI_MAX, 0, comm);
        t_bal_g[1] = t_bal_g[1] / npes;

        localSz=pNodesBalanced.size();
        par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);
        numBalOct=globalSz;

        if(!rank) {
            std::cout << YLW<<"Number of balanced octants: " << numBalOct <<NRM<< std::endl;
            std::cout<< YLW<<" Tree balancing time (max): "<<t_bal_g[2]<<NRM<<std::endl;
        }
        MPI_Barrier(comm);

        if (!rank) std::cout << RED << "mesh generation begin" << NRM << std::endl;
        t1 = MPI_Wtime();//std::chrono::high_resolution_clock::now();
        ot::Mesh mesh(pNodesBalanced, stencilSz, eleOrder,comm,false,ot::FEM_CG);
        t2 = MPI_Wtime();//std::chrono::high_resolution_clock::now();
        t_mesh = t2-t1;//std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        if (!rank) std::cout << RED << "mesh generation end" << NRM << std::endl;

        par::Mpi_Reduce(&t_mesh, t_mesh_g, 1, MPI_MIN, 0, comm);
        par::Mpi_Reduce(&t_mesh, t_mesh_g + 1, 1, MPI_SUM, 0, comm);
        par::Mpi_Reduce(&t_mesh, t_mesh_g + 2, 1, MPI_MAX, 0, comm);
        t_mesh_g[1] = t_mesh_g[1] / npes;

        if(!rank) {
           std::cout<< YLW<<"Mesh generation time (max): "<<t_mesh_g[2]<<NRM<<std::endl;
                std::cout<<"\t"<<YLW<<" e2e (min,mean,max): "<<"( "<<t_e2e_g[0]<<"\t"<<t_e2e_g[1]<<"\t"<<t_e2e_g[2]<<" )"<<NRM<<std::endl;
                std::cout<<"\t"<<YLW<<" e2n (min,mean,max): "<<"( "<<t_e2n_g[0]<<"\t"<<t_e2n_g[1]<<"\t"<<t_e2n_g[2]<<" )"<<NRM<<std::endl;
                std::cout<<"\t"<<YLW<<" sm (min,mean,max): "<<"( "<<t_sm_g[0]<<"\t"<<t_sm_g[1]<<"\t"<<t_sm_g[2]<<" )"<<NRM<<std::endl;
                std::cout<<"\t"<<YLW<<" blk (min,mean,max): "<<"( "<<t_blk_g[0]<<"\t"<<t_blk_g[1]<<"\t"<<t_blk_g[2]<<" )"<<NRM<<std::endl;
        }

        localSz = mesh.getDegOfFreedom();
        par::Mpi_Reduce(&localSz,&cg_sz,1,MPI_SUM,0,comm);


    }else
    {
        // sequential run:
        //SFC::seqSort::SFC_treeSort(&(*tmpNodes.begin()),tmpNodes.size(),pNodesSorted,pNodesConstructed,pNodesBalanced,maxDepth,maxDepth,root,ROOT_ROTATION,1,TS_CONSTRUCT_OCTREE);

        if (!rank) std::cout << RED << "Remove duplicates begin" << NRM << std::endl;
        auto t1=MPI_Wtime();
        SFC::seqSort::SFC_treeSort(&(*(tmpNodes.begin())),tmpNodes.size(),pNodesSorted,pNodesConstructed,pNodesBalanced,maxDepth,maxDepth,root,ROOT_ROT_ID,1,TS_REMOVE_DUPLICATES);
        auto t2=MPI_Wtime();
        t_rd=t2-t1;
        if (!rank) std::cout << RED << "Remove duplicates end" << NRM << std::endl;
        tmpNodes.clear();

        localSz=pNodesSorted.size();
        par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);
        numUniqueOct=globalSz;

        t_rd_g[0]=t_rd;
        t_rd_g[1]=t_rd;
        t_rd_g[2]=t_rd;

        if(!rank)
        {
            std::cout << YLW<< "Number of unique octants: "<<numUniqueOct<<NRM<<std::endl;
            std::cout << YLW<< "Remove duplicates time(max): "<<t_rd_g[2]<<NRM<<std::endl;
        }

        if (!rank) std::cout << RED << "Tree construction begin" << NRM << std::endl;
        t1=MPI_Wtime();
        SFC::seqSort::SFC_treeSort(&(*(pNodesSorted.begin())),pNodesSorted.size(),pNodesSorted,pNodesConstructed,pNodesBalanced,maxDepth,maxDepth,root,ROOT_ROT_ID,1,TS_CONSTRUCT_OCTREE);
        t2=MPI_Wtime();
        t_cons=t2-t1;
        if (!rank) std::cout << RED << "Tree construction end" << NRM << std::endl;
        pNodesSorted.clear();

        t_cons_g[0]=t_cons;
        t_cons_g[1]=t_cons;
        t_cons_g[2]=t_cons;


        localSz=pNodesConstructed.size();
        par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);
        numConsOct=globalSz;
        if(!rank)
        {
            std::cout << YLW<< "Number of unbalanced octants: "<<numConsOct<<NRM<<std::endl;
            std::cout << YLW<< "Tree construction time(max): "<<t_cons_g[2]<<NRM<<std::endl;
        }


        if (!rank) std::cout << RED << "2:1 balance begin" << NRM << std::endl;
        t1 = MPI_Wtime();//std::chrono::high_resolution_clock::now();
        SFC::seqSort::SFC_treeSort(&(*(pNodesConstructed.begin())),pNodesConstructed.size(),pNodesSorted,pNodesConstructed,pNodesBalanced,maxDepth,maxDepth,root,ROOT_ROT_ID,1,TS_BALANCE_OCTREE);
        t2 = MPI_Wtime();//std::chrono::high_resolution_clock::now();
        t_bal = t2-t1;//std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        if (!rank) std::cout << RED << "2:1 balance end" << NRM << std::endl;

        t_bal_g[0]=t_bal;
        t_bal_g[1]=t_bal;
        t_bal_g[2]=t_bal;

        localSz=pNodesBalanced.size();
        par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);
        numBalOct=globalSz;

        if(!rank) {
            std::cout << YLW<<"Number of balanced octants: " << numBalOct <<NRM<< std::endl;
            std::cout<< YLW<<" Tree balancing time (max): "<<t_bal_g[2]<<NRM<<std::endl;
        }


        // Note that is balance test works only if you have the global tree.
/*        if(!ot::test::isBalanced(m_uiDim,m_uiMaxDepth,"failedOct",pNodesBalanced,true,1))
        {
           std::cout<<"Balance test failed. "<<std::endl;
        }*/

        if (!rank) std::cout << RED << "mesh generation begin" << NRM << std::endl;
        t1 = MPI_Wtime();//std::chrono::high_resolution_clock::now();
        {
            ot::Mesh mesh(pNodesBalanced,stencilSz,eleOrder,comm,false,ot::FEM_CG);
            cg_sz = mesh.getDegOfFreedom();
        }
        t2 = MPI_Wtime();//std::chrono::high_resolution_clock::now();
        t_mesh = t2-t1;//std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        if (!rank) std::cout << RED << "mesh generation end" << NRM << std::endl;

        t_mesh_g[0]=t_mesh;
        t_mesh_g[1]=t_mesh;
        t_mesh_g[2]=t_mesh;



        if(!rank) {
            std::cout<< YLW<<"Mesh generation time (max): "<<t_mesh_g[2]<<NRM<<std::endl;
            std::cout<<"\t"<<YLW<<" e2e (min,mean,max): "<<"( "<<t_e2e_g[0]<<"\t"<<t_e2e_g[1]<<"\t"<<t_e2e_g[2]<<" )"<<NRM<<std::endl;
            std::cout<<"\t"<<YLW<<" e2n (min,mean,max): "<<"( "<<t_e2n_g[0]<<"\t"<<t_e2n_g[1]<<"\t"<<t_e2n_g[2]<<" )"<<NRM<<std::endl;
            std::cout<<"\t"<<YLW<<" sm (min,mean,max): "<<"( "<<t_sm_g[0]<<"\t"<<t_sm_g[1]<<"\t"<<t_sm_g[2]<<" )"<<NRM<<std::endl;
            std::cout<<"\t"<<YLW<<" blk (min,mean,max): "<<"( "<<t_blk_g[0]<<"\t"<<t_blk_g[1]<<"\t"<<t_blk_g[2]<<" )"<<NRM<<std::endl;
        }


        dollar::text(std::cout);

        


    }
    if(!rank)
    {
        char fileName[256];
        sprintf(fileName,"%s_%d",prefix,npes);
        //std::cout<<"filename: "<<fileName<<std::endl;
        std::ofstream statFile;
        statFile.open(fileName);

        if(npes==2)
        {
            statFile<<"npes\tgrainSz\tdim\tmaxDepth\ttolerance\tsf_k\tstencilSz\teleOrder\tt_rm_min\tt_rm_mean\tt_rm_max\tt_cons_min\tt_cons_mean\tt_cons_max\tt_bal_min\tt_bal_mean\tt_bal_max\tt_mesh_min\tt_mesh_mean\tt_mesh_max\tnum_oct\tnum_oct_cons\tnum_oct_bal\tcg_nodes"<<std::endl;
        }

        statFile<<npes<<"\t"<<numPts<<"\t"<<dim<<"\t"<<maxDepth<<"\t"<<tol<<"\t"<<sf_k<<"\t"<<stencilSz<<"\t"<<eleOrder<<"\t"<<t_rd_g[0]<<"\t"<<t_rd_g[1]<<"\t"<<t_rd_g[2]<<"\t"<<t_cons_g[0]<<"\t"<<t_cons_g[1]<<"\t"<<t_cons_g[2]<<"\t"<<t_bal_g[0]<<"\t"<<t_bal_g[1]<<"\t"<<t_bal_g[2]<<"\t"<<t_mesh_g[0]<<"\t"<<t_mesh_g[1]<<"\t"<<t_mesh_g[2]<<"\t"<<numUniqueOct<<"\t"<<numConsOct<<"\t"<<numBalOct<<"\t"<<cg_sz<<std::endl;
        statFile.close();
    }

}


void weakScalingDriver(char * ptsFile,bool genPts,DendroIntL numPts,unsigned int dim,unsigned int maxDepth,double tolerance,int distribution,unsigned int k,unsigned int options,unsigned int sf_k,char * prefix,MPI_Comm comm)
{

    int rank,npes;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    if(!rank) std::cout<<"======================================================================"<<std::endl;
    if(!rank) std::cout<<"     Warm Up run  (par) begin                               "<<std::endl;
    if(!rank) std::cout<<"======================================================================"<<std::endl;

    meshBenchMark(ptsFile,genPts,numPts,dim,maxDepth,distribution,tolerance,sf_k,options,prefix,comm);

    if(!rank) std::cout<<"======================================================================"<<std::endl;
    if(!rank) std::cout<<"     Warm Up run  (par) end                               "<<std::endl;
    if(!rank) std::cout<<"======================================================================"<<std::endl;

    MPI_Barrier(comm);

    // full run
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    if(!rank) std::cout<<"======================================================================"<<std::endl;
    if(!rank) std::cout<<"     Beginning Full Run.   "<<std::endl;
    if(!rank) std::cout<<"======================================================================"<<std::endl;

    meshBenchMark(ptsFile,genPts,numPts,dim,maxDepth,distribution,tolerance,sf_k,options,prefix,comm);

#ifdef RUN_WEAK_SCALING

    int proc_group = 0;
    int min_np = 2;
    for (int i = npes; rank < i && i >= min_np; i = i >> 1) proc_group++;
    MPI_Comm comm_ws;

    MPI_Comm_split(comm, proc_group, rank, &comm_ws);

    MPI_Comm_rank(comm_ws, &rank);
    MPI_Comm_size(comm_ws, &npes);

    MPI_Barrier(comm_ws);
    meshBenchMark(ptsFile,genPts,numPts,dim,maxDepth,distribution,tolerance,sf_k,options,prefix,comm_ws);


    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);
    //MPI_Barrier(comm);

    if(!rank) std::cout<<"======================================================================"<<std::endl;
    if(!rank) std::cout<<"     Weak Scaling Run Complete.   "<<std::endl;
    if(!rank) std::cout<<"======================================================================"<<std::endl;

#endif




}



int main(int argc, char** argv) {

    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);

    int rank, npes;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);

    unsigned int options = TS_BALANCE_OCTREE;


    if (argc < 3) {
        if (!rank)
            std::cout << "Usage :" << argv[0] << " inp  numPts " << " dim " << " maxDepth " << " tol  "<< " distribution (default =0) (0- Normal 1- Uniform 2- LogarithmicNormal) genPts(default=1) SplitterFixK (default =2)"
                      << std::endl;
        return -1;
    }

    DendroIntL numPts = atoll(argv[2]);
    unsigned int dim = atoi(argv[3]);
    unsigned int maxDepth = atoi(argv[4]);
    m_uiMaxDepth=maxDepth;
    double tol = 0.001;
    unsigned int distribution = 0;
    unsigned int sf_k = 2;

    unsigned int k_perOct=1;
    bool genPts=true;

    std::string prefix="meshProf";


    if (argc > 5)
        tol = atof(argv[5]);

    if (argc > 6)
        distribution = atoi(argv[6]);

    if (argc > 7)
         genPts=(bool)atoi(argv[7]);

    if (argc > 8)
        sf_k = atoi(argv[8]);

    if (!rank) std::cout << "sf parameter: " << sf_k << std::endl;

    _InitializeHcurve(dim);
    if (!rank) std::cout << "Initialized H-Curves for dimension " << dim << std::endl;

    weakScalingDriver(argv[1],genPts,numPts,dim,maxDepth,tol,distribution,k_perOct,options,sf_k,(char* )prefix.c_str(),MPI_COMM_WORLD);
    dollar::text(std::cout);
    MPI_Finalize();

    return 0;

}

