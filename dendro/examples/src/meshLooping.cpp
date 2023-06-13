//
// Created by milinda on 12/15/16.
//

/**
 * @author: Milinda Fernando
 * @breif: Contains some example looping over elelements in wavelet DA
 *
 * */


#include "TreeNode.h"
#include "mpi.h"
#include "genPts_par.h"
#include "sfcSort.h"
#include "mesh.h"
#include "dendro.h"
#define ROOT_ROT_ID 0

int main(int argc, char** argv) {

    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);
    int rank, npes;
    MPI_Comm comm=MPI_COMM_WORLD;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);



    if (argc < 3) {
        if (!rank)
            std::cout << "Usage :" << argv[0] << " numPts " << " dim " << " maxDepth " << " tol  "
                      << " distribution (default =0) (0- Normal 1- Uniform 2- LogarithmicNormal) SplitterFixK (default =2)"
                      << std::endl;
    }

    DendroIntL numPts = atoll(argv[1]);
    unsigned int dim = atoi(argv[2]);
    unsigned int maxDepth = atoi(argv[3]);
    m_uiMaxDepth=maxDepth;
    double tol = 0.001;
    unsigned int distribution = 0;
    unsigned int sf_k = 2;

    unsigned int k_perOct=1;


    std::vector<ot::TreeNode> pNodesSorted;
    std::vector<ot::TreeNode> pNodesConstructed;
    std::vector<ot::TreeNode> pNodesBalanced;

    std::vector<ot::TreeNode> tmpNodes;
    m_uiMaxDepth=maxDepth; // this is needed because we increase this when we embedded the tree.
    unsigned int stencilSz=1;

    DendroIntL totPts=numPts*dim;
    double * pts=new double[totPts];


    if (argc > 4)
        tol = atof(argv[4]);

    if (argc > 5)
        distribution = atoi(argv[5]);

    if (argc > 6)
        sf_k = atoi(argv[6]);

    if (!rank) std::cout << "sf parameter: " << sf_k << std::endl;

    _InitializeHcurve(dim);
    if (!rank) std::cout << "Initialized H-Curves for dimension " << dim << std::endl;


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

    if(distribution==0)
        genGauss(0.15,numPts,dim,pts);
    else if(distribution==1)
        genUniformRealDis(0.1,numPts,dim,pts);
    else if(distribution==2)
        genLogarithmicGauss(0.5,numPts,dim,pts);
    else
        genGauss(0.5,numPts,dim,pts);  // Default case.


    pts2Octants(tmpNodes,pts,totPts,dim,maxDepth); // generate octants from the points.

    pNodesSorted.clear();
    pNodesConstructed.clear();
    pNodesBalanced.clear();
    ot::TreeNode root(0,0,0,0,dim,maxDepth);
    double t_g=0;
    double t_bal=0;
    double t_mesh=0;

    if(npes>=2) {

        if (!rank) std::cout << RED << "2:1 balance begin" << NRM << std::endl;
        auto t1 = MPI_Wtime();//std::chrono::high_resolution_clock::now();
        SFC::parSort::SFC_treeSort(tmpNodes, pNodesSorted, pNodesConstructed, pNodesBalanced, tol, maxDepth, root,
                                   ROOT_ROT_ID, 1, TS_BALANCE_OCTREE, sf_k, comm);
        auto t2 = MPI_Wtime();//std::chrono::high_resolution_clock::now();
        t_bal = t2-t1;//std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        if (!rank) std::cout << RED << "2:1 balance end" << NRM << std::endl;
        par::Mpi_Reduce(&t_bal,&t_g,1,MPI_MAX,0,comm);
        DendroIntL localSz=pNodesBalanced.size();
        DendroIntL globalSz;
        par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);

        if(!rank)
        {
            std::cout<<"octree construction + balance  time: "<<t_g<<std::endl;
            std::cout<<"Number of bal octs: "<<globalSz<<std::endl;
        }



        if(!rank)
            std::cout<<"Number of balanced octants: "<<globalSz<<std::endl;
        assert(par::test::isUniqueAndSorted(pNodesBalanced,comm));

        if (!rank) std::cout << RED << "mesh generation begin" << NRM << std::endl;
        t1 = MPI_Wtime();//std::chrono::high_resolution_clock::now();
        ot::Mesh mesh(pNodesBalanced, stencilSz,1, comm);
        t2 = MPI_Wtime();//std::chrono::high_resolution_clock::now();
        t_mesh = t2-t1;//std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        if (!rank) std::cout << RED << "mesh generation end" << NRM << std::endl;
        par::Mpi_Reduce(&t_mesh,&t_g,1,MPI_MAX,0,comm);

        if(!rank)
          std::cout<<"Mesh time: "<<t_g<<std::endl;

      /*  for(mesh.init<ot::WaveletDA::INDEPENDENT>();mesh.nextAvailable<ot::WaveletDA::INDEPENDENT>();mesh.next<ot::WaveletDA::INDEPENDENT>())
        {

        }*/


    }else
    {
        if (!rank) std::cout << RED << "2:1 balance begin" << NRM << std::endl;
        auto t1 = MPI_Wtime();//std::chrono::high_resolution_clock::now();
        SFC::parSort::SFC_treeSort(tmpNodes, pNodesSorted, pNodesConstructed, pNodesBalanced, tol, maxDepth, root,
                                   ROOT_ROT_ID, 1, TS_BALANCE_OCTREE, sf_k, comm);
        auto t2 = MPI_Wtime();//std::chrono::high_resolution_clock::now();
        t_bal = t2-t1;//std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        if (!rank) std::cout << RED << "2:1 balance end" << NRM << std::endl;
        par::Mpi_Reduce(&t_bal,&t_g,1,MPI_MAX,0,comm);
        DendroIntL localSz=pNodesBalanced.size();
        DendroIntL globalSz;
        par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);

        if(!rank)
        {
            std::cout<<"octree construction + balance  time: "<<t_g<<std::endl;
            std::cout<<"Number of bal octs: "<<globalSz<<std::endl;
        }



        if(!rank)
            std::cout<<"Number of balanced octants: "<<globalSz<<std::endl;

        if (!rank) std::cout << RED << "mesh generation begin" << NRM << std::endl;
        t1 = MPI_Wtime();//std::chrono::high_resolution_clock::now();
        ot::Mesh mesh(pNodesBalanced, stencilSz,1,comm);
        t2 = MPI_Wtime();//std::chrono::high_resolution_clock::now();
        t_mesh = t2-t1;//std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        if (!rank) std::cout << RED << "mesh generation end" << NRM << std::endl;
        par::Mpi_Reduce(&t_mesh,&t_g,1,MPI_MAX,0,comm);

        if(!rank)
            std::cout<<"Mesh time: "<<t_g<<std::endl;
    }



    MPI_Finalize();

    return 0;

}