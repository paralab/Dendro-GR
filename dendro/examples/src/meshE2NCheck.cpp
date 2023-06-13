//
// Created by milinda on 2/25/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief The purpose of this file is to check the E2N mapping , is correct. We generate a complete octree and compute number of tasks using p number of processors, then number of actual variables should be
 * equal to the number of actual variables in the sequential case.
*/
//
#define ROOT_ROT_ID 0



#include "TreeNode.h"
#include "mpi.h"
#include "genPts_par.h"
#include "sfcSort.h"
#include "mesh.h"
#include "dendro.h"
#include "dendroIO.h"
#include "dollar.hpp"
#include "octUtils.h"


int main (int argc, char** argv)
{


    MPI_Init(&argc,&argv);

    MPI_Comm comm=MPI_COMM_WORLD;

    int rank,npes;
    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&npes);

    unsigned int options = TS_BALANCE_OCTREE;


    if (argc < 3) {
        if (!rank)
            std::cout << "Usage :" << argv[0] << " inp  numPts " << " dim " << " maxDepth " << " tol  "<< " distribution (default =0) (0- Normal 1- Uniform 2- LogarithmicNormal) genPts(default=1) SplitterFixK (default =2)"
                      << std::endl;
        return -1;
    }
    char * ptsFile=argv[1];
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

    if(npes==1)
    {
        DendroIntL totPts=numPts*dim;
        char pFile[256];
        std::vector<double> ptsFromFile;
        sprintf(pFile, "%s%d_%d.pts", ptsFile, rank, npes);
        double * pts=new double[totPts];
        if(genPts) genGauss(0.15,numPts,dim,ptsFile,comm);

        std::vector<ot::TreeNode> tmpNodes;
        m_uiMaxDepth=maxDepth;
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

        delete [] pts;
        /*std::set<ot::Key,OctreeComp<ot::Key>> tmpKeySet;

        std::vector<ot::Key> tmpKeys;
        std::vector<ot::Key> tmpKeys1;
        ot::Key rootKey(0,0,0,0,m_uiDim,m_uiMaxDepth);

        for(unsigned int i=0;i<tmpNodes.size();i++)
        {
           *//* tmpKeys.push_back(ot::Key(tmpNodes[i]));
            tmpKeys.back().addOwner(i);
            if(i%2==0)
            {
                tmpKeys.back().addOwner(i+1);
            }*//*

            auto it=tmpKeySet.emplace(ot::Key(tmpNodes[i]));
            it.first->addOwner(i);
            if(i%2==0)
                it.first->addOwner(i+1);
            if(it.second)
                tmpKeys.push_back(*(it.first));

        }

        //tmpKeys.insert(tmpKeys.end(),tmpKeySet.begin(),tmpKeySet.end());
        std::cout<<"begin SFC Sort "<<tmpKeys.size()<<std::endl;
        SFC::seqSort::SFC_treeSort(&(*(tmpKeys.begin())),tmpKeys.size(),tmpKeys1,tmpKeys1,tmpKeys1,m_uiMaxDepth,m_uiMaxDepth,rootKey,0,1,TS_SORT_ONLY);
        assert(seq::test::isSorted(tmpKeys));
        ot::Key minKey=tmpKeys.front();
        ot::Key minKey1;
        ot::Key maxKey1;
        ot::Key maxKey=tmpKeys.back();
        SFC::seqSort::SFC_treeSortLocalOptimal(&(*(tmpKeys.begin())),tmpKeys.size(),m_uiMaxDepth,m_uiMaxDepth,rootKey,ROOT_ROTATION,true,minKey1);
        SFC::seqSort::SFC_treeSortLocalOptimal(&(*(tmpKeys.begin())),tmpKeys.size(),m_uiMaxDepth,m_uiMaxDepth,rootKey,ROOT_ROTATION,false,maxKey1);
        if(minKey!=minKey1)
        {
            std::cout<<"minKey: "<<minKey<<" minKey1: "<<minKey1<<std::endl;
        }
        assert(minKey==minKey1);
        assert(maxKey==maxKey1);


        std::cout<<"tmpKeys Sorted Size: "<<tmpKeys.size()<<std::endl;*/
        if(!rank) std::cout<<"This test is intended to run in parallel. "<<std::endl;

       /* RefElement m_uiRefEl;

        m_uiRefEl=RefElement(m_uiDim,1);
        m_uiRefEl.generateHeaderFile("order_1.h");
        m_uiRefEl=RefElement(m_uiDim,2);
        m_uiRefEl.generateHeaderFile("order_2.h");
        m_uiRefEl=RefElement(m_uiDim,4);
        m_uiRefEl.generateHeaderFile("order_4.h");*/

        return 0;
    }else
    {

        std::vector<ot::TreeNode> pNodesSorted;
        std::vector<ot::TreeNode> pNodesConstructed;
        std::vector<ot::TreeNode> pNodesBalanced;
        std::vector<ot::TreeNode> pNodesBalanced_global;

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




        unsigned int stencilSz=1;

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
            par::Mpi_Allreduce(&localSz,&globalSz,1,MPI_SUM,comm);
            numBalOct=globalSz;

            if(!rank) {
                std::cout << YLW<<"Number of balanced octants: " << numBalOct <<NRM<< std::endl;
                std::cout<< YLW<<" Tree balancing time (max): "<<t_bal_g[2]<<NRM<<std::endl;
            }
            MPI_Barrier(comm);


          /*  std::vector<ot::Block> blockList;
            unsigned int d_min,d_max;
            octree2BlockDecomposition(pNodesBalanced,blockList,m_uiMaxDepth,d_min,d_max,0,pNodesBalanced.size());
            //std::cout<<"rank: "<<rank<<" block list size: "<<blockList.size()<<" d_min: "<<d_min<<" d_max: "<<d_max<<std::endl;

            assert(ot::test::isBlockListValid(pNodesBalanced,blockList,d_min,d_max,0,pNodesBalanced.size()));*/

            //Need to gather the balanced octree to perform the test.
            /*unsigned int localBalSz=(pNodesBalanced.size())*sizeof(ot::TreeNode);
            unsigned int * recvCount=new unsigned int [npes];
            unsigned int * recvCountOffset=new unsigned int [npes];

            par::Mpi_Allgather(&localBalSz,recvCount,1,comm);
            recvCountOffset[0]=0;
            omp_par::scan(recvCount,recvCountOffset,npes);

            pNodesBalanced_global.resize((recvCount[npes-1]+recvCountOffset[npes-1])/sizeof(ot::TreeNode));


            MPI_Allgatherv(&(*(pNodesBalanced.begin())),localBalSz,MPI_BYTE,&(*(pNodesBalanced_global.begin())),(int *)recvCount,(int *)recvCountOffset,MPI_BYTE,comm);

            delete [] recvCount;
            delete [] recvCountOffset;
            //if(!rank) treeNodesTovtk(pNodesBalanced_global,rank,"balOct_g");
            assert(pNodesBalanced_global.size()==numBalOct);*/

            if (!rank) std::cout << RED << "mesh generation begin" << NRM << std::endl;
            t1 = MPI_Wtime();//std::chrono::high_resolution_clock::now();
            ot::Mesh mesh(pNodesBalanced, stencilSz, 4,comm);
            t2 = MPI_Wtime();//std::chrono::high_resolution_clock::now();
            t_mesh = t2-t1;//std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
            if (!rank) std::cout << RED << "mesh generation end" << NRM << std::endl;

            par::Mpi_Reduce(&t_mesh, t_mesh_g, 1, MPI_MIN, 0, comm);
            par::Mpi_Reduce(&t_mesh, t_mesh_g + 1, 1, MPI_SUM, 0, comm);
            par::Mpi_Reduce(&t_mesh, t_mesh_g + 2, 1, MPI_MAX, 0, comm);
            t_mesh_g[1] = t_mesh_g[1] / npes;
            unsigned int localNodes=mesh.getNumLocalMeshNodes();
            unsigned int localNodes_g=0;
            par::Mpi_Allreduce(&localNodes,&localNodes_g,1,MPI_SUM,comm);


            if(!rank) {
                std::cout<< YLW<<"Mesh generation time (max): "<<t_mesh_g[2]<<NRM<<std::endl;
            }

#ifdef SEQUENTIAL_OCTREE_MESH_COMPARE
                m_uiMaxDepth=maxDepth; // need to reset the maxdepth for the octree.
                if(!rank)std::cout<<"Staring sequential run to see whether we get the same number of local nodes "<<std::endl;
                ot::Mesh mesh1(pNodesBalanced_global, stencilSz, 1);

                mesh.init<ot::WaveletDA::INDEPENDENT>();
                mesh1.init<ot::WaveletDA::INDEPENDENT>();
                unsigned int* eleNeig=new unsigned int[mesh.getNumDirections()];
                unsigned int* eleNeig1=new unsigned int[mesh.getNumDirections()];

                unsigned int* nodeNeigh=new unsigned int [mesh.getNumNodesPerElement()];
                unsigned int* nodeNeigh1=new unsigned int [mesh.getNumNodesPerElement()];


                const ot::TreeNode* elePtr=mesh.getPtrAllElements();
                const ot::TreeNode* elePtr1=mesh1.getPtrAllElements();
                unsigned int x,y,z,sz,x1,y1,z1,sz1;
                unsigned int e,i,j,k;
                unsigned int e1,i1,j1,k1;

                std::vector<ot::TreeNode> currentEle;
                std::vector<ot::TreeNode> owner;
                std::vector<ot::TreeNode> owner1;




                while(mesh1.nextAvailable<ot::WaveletDA::INDEPENDENT>() &&(elePtr[mesh.currentIndex()]!=elePtr1[mesh1.currentIndex()])) {mesh1.next<ot::WaveletDA::INDEPENDENT>();}

                    while (mesh.nextAvailable<ot::WaveletDA::INDEPENDENT>() &&
                           mesh1.nextAvailable<ot::WaveletDA::INDEPENDENT>()) {

                        //std::cout<<"rank : "<<rank<<" checking: "<<mesh.currentIndex()<<" and "<<mesh1.currentIndex()<<std::endl;
                        assert(elePtr[mesh.currentIndex()]==elePtr1[mesh1.currentIndex()]);
                        mesh.currentElementNeighbourIndexList(eleNeig);
                        mesh1.currentElementNeighbourIndexList(eleNeig1);

                        mesh.currentElementNodeList_DG(nodeNeigh);
                        mesh1.currentElementNodeList_DG(nodeNeigh1);

                        for (unsigned int dir = 0; dir < mesh.getNumDirections(); dir++) {
                            if (eleNeig[dir] == UINT_MAX) {
                                assert(eleNeig1[dir] == UINT_MAX);
                            } else {
                                if(elePtr[eleNeig[dir]] != elePtr1[eleNeig1[dir]])
                                {
                                    std::cout<<"rank: "<<rank<<" E2E mismatch: "<<dir<<" element: "<<elePtr[mesh.currentIndex()]<<" & element: "<<elePtr1[mesh1.currentIndex()]<<" neighbour: "<<elePtr[eleNeig[dir]]<<" & "<<elePtr1[eleNeig1[dir]]<<std::endl;
                                }
                                assert(elePtr[eleNeig[dir]] == elePtr1[eleNeig1[dir]]);
                            }
                        }


                         // NOTE: !!! To use the following test E2N mapping should be in DG indexing after removing duplicates
                           /* for(unsigned int node=0;node<mesh.getNumNodesPerElement();node++)
                            {


                                mesh.dg2eijk(nodeNeigh[node],e,i,j,k);
                                mesh1.dg2eijk(nodeNeigh1[node],e1,i1,j1,k1);

                                x=elePtr[e].getX();
                                y=elePtr[e].getY();
                                z=elePtr[e].getZ();
                                sz=1u<<(m_uiMaxDepth-elePtr[e].getLevel());

                                x1=elePtr1[e1].getX();
                                y1=elePtr1[e1].getY();
                                z1=elePtr1[e1].getZ();
                                sz1=1u<<(m_uiMaxDepth-elePtr1[e1].getLevel());

                                if(((x1+i1*sz1)!=(x+i*sz)) || ((y1+j1*sz1)!=(y+j*sz))|| ((z1+k1*sz1)!=(z+k*sz))){
                                    currentEle.push_back(elePtr[mesh.currentIndex()]);
                                    owner.push_back(elePtr[e]);
                                    owner1.push_back(elePtr1[e1]);
                                    std::cout<<"[E2N Mapping Failed: ] rank "<<rank <<"E2N failed for :" <<elePtr[mesh.currentIndex()]<<" node : "<<node<<" element e1: "<<e1<<" : "<<elePtr1[e1] <<" element e: "<<e<<" : "<<elePtr[e]<<" x1: "<<(x1+i1*sz1)<<" x: "<<(x+i*sz)<< " y1: "<<(y1+i1*sz1)<<" y: "<<(y+i*sz)<< " z1: "<<(z1+i1*sz1)<<" z: "<<(z+i*sz)<<std::endl;

                                    std::cout<<"E2E mapping for failed element: "<<elePtr[mesh.currentIndex()]<<": ";
                                    for (unsigned int dir = 0; dir < mesh.getNumDirections(); dir++)
                                    {
                                        std::cout<<" ("<<elePtr[eleNeig[dir]]<<", "<<elePtr1[eleNeig1[dir]]<<" ) ";

                                    }
                                    std::cout<<std::endl;

                                }



                            }*/

                        mesh.next<ot::WaveletDA::INDEPENDENT>();
                        mesh1.next<ot::WaveletDA::INDEPENDENT>();


                    }

                // to write the invalid E2N mapping entries.
                /*treeNodesTovtk(currentEle,rank,"currentElement");
                treeNodesTovtk(owner,rank,"owner");
                treeNodesTovtk(owner1,rank,"owner1");*/

                delete[] eleNeig;
                delete[] eleNeig1;
                delete[] nodeNeigh;
                delete[] nodeNeigh1;

                // summation of number of local nodes in each processor should be equal to the number of local nodes in the sequential case.
                if(!rank) std::cout<<"seq: mesh local Nodes: "<<mesh1.getNumLocalMeshNodes()<<" par: mesh local nodes: "<<localNodes_g<<std::endl;
                assert(mesh1.getNumLocalMeshNodes()==localNodes_g);

#endif



    }

    MPI_Finalize();





    return 0;
}

