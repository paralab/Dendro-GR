//
// Created by milinda on 8/30/18.
//

/**
 * @author Milinda Fernando,
 * School of Computing, University of Utah
 * @brief generate an octree balance it, and perform comparison study between binary search and sfc search.
 *
 **/


#include "TreeNode.h"
#include "sfcSort.h"
#include "sfcSearch.h"
#include "genPts_par.h"
#include "key.h"
#include "skey.h"


#include <iostream>
#include <vector>
#include <skey.h>





int main(int argc, char** argv) {

    MPI_Init(&argc, &argv);

    int rank, npes;
    MPI_Comm GLOBAL_COMM = MPI_COMM_WORLD;
    MPI_Comm_rank(GLOBAL_COMM, &rank);
    MPI_Comm_size(GLOBAL_COMM, &npes);

    if (argc < 3) {
        if (!rank)
            std::cout << "Usage :" << argv[0] << " numPts numKeys" << " dim " << " maxDepth " << " tol[0.1]  "
                      << " distribution (0- Normal 1- Uniform 2- LogarithmicNormal)" << std::endl;
    }

    DendroIntL numPts = atoll(argv[1]);
    DendroIntL numKeys = atoll(argv[2]);
    unsigned int dim = atoi(argv[3]);
    unsigned int maxDepth = atoi(argv[4]);
    m_uiMaxDepth=maxDepth;
    double tol = 0.1;
    unsigned int distribution = 0;
    if (argc > 5)
        tol = atof(argv[5]);

    if (argc > 6)
        distribution = atoi(argv[6]);


    _InitializeHcurve(dim);
    //if (!rank) std::cout << "Initialized Hcurves for dimention " << dim << std::endl;

    DendroIntL localSz, totalSz;

    DendroIntL totPts=numPts*dim;
    double * pts=new double[totPts];

    if(distribution==0)
        genGauss(0.1,numPts,dim,pts);
    else if(distribution==1)
        genUniformRealDis(0.1,numPts,dim,pts);
    else if(distribution==2)
        genLogarithmicGauss(0.5,numPts,dim,pts);
    else
        genGauss(0.5,numPts,dim,pts);  // Default case.

    std::vector<ot::TreeNode> tmpNodes;

    pts2Octants(tmpNodes,pts,totPts,dim,maxDepth); // generate octants from the points.

    delete [] pts;

    //std::cout<<"inp size: "<<tmpNodes.size()<<std::endl;

    std::vector<ot::TreeNode> pSorted;
    std::vector<ot::TreeNode> pConstructed;
    std::vector<ot::TreeNode> pBalanced;

    ot::TreeNode root(0,0,0,0,dim,maxDepth);


    localSz=tmpNodes.size();
    par::Mpi_Reduce(&localSz,&totalSz,1,MPI_SUM,0,GLOBAL_COMM);

    if(!rank) std::cout<<YLW<<"================================================================="<<NRM<<std::endl;
    if(!rank) std::cout<<YLW<<" Total number of Points: "<<totalSz<<NRM<<std::endl;
    if(!rank) std::cout<<YLW<<"================================================================="<<NRM<<std::endl;

    double t1, t2;


    if(npes>2)
    {
        std::cout<<"This is a sequential test"<<std::endl;
        return 0;
    }


    if(npes==1)
    {
        if(!rank) std::cout <<RED<<"==================================================================================================="<<NRM<<std::endl;
        if(!rank) std::cout <<RED<<"           Executing the SEQUENTIAL TREESORT with options : "<<NRM<<std::endl;
        if(!rank) std::cout<<RED<<"==================================================================================================="<<NRM<<std::endl;


        t1=MPI_Wtime();
        SFC::seqSort::SFC_treeSort(&(*(tmpNodes.begin())),tmpNodes.size(),pSorted,pConstructed,pBalanced,maxDepth,maxDepth,root,0,1,TS_CONSTRUCT_OCTREE);
        t2=MPI_Wtime();

        if(!rank) std::cout<<"construction time(s) : "<<(t2-t1)<<std::endl;
        if(!rank) std::cout<<"number of octs(construction): "<<pConstructed.size()<<std::endl;

        t1=MPI_Wtime();
        SFC::seqSort::SFC_treeSort(&(*(pConstructed.begin())),pConstructed.size(),pSorted,pConstructed,pBalanced,maxDepth,maxDepth,root,0,1,TS_BALANCE_OCTREE);
        t2=MPI_Wtime();

        if(!rank) std::cout<<"balancing time(s) : "<<(t2-t1)<<std::endl;
        if(!rank) std::cout<<"number of octs(balOcts): "<<pBalanced.size()<<std::endl;

        assert(seq::test::isUniqueAndSorted(pSorted));
        assert(seq::test::isUniqueAndSorted(pConstructed));
        assert(seq::test::isUniqueAndSorted(pBalanced));

        pSorted.clear();
        pConstructed.clear();

        std::vector<ot::TreeNode> tmpKeys;
        std::vector<unsigned int > tmpKeyResults;
        tmpNodes.clear();

        const unsigned int NUM_ITERARTIONS=10;

        double * pts_keys=new double[numKeys*dim];


        // warm up run.

        tmpNodes.clear();
        tmpKeys.clear();
        tmpKeyResults.clear();
        genGauss(0.1,numKeys,dim,pts_keys);
        pts2Octants(tmpNodes,pts_keys,numKeys*dim,dim,m_uiMaxDepth);

        for(unsigned int i=0;i<tmpNodes.size();i++)
            tmpKeys.push_back(ot::TreeNode(tmpNodes[i]));


        //SFC::seqSearch::SFC_treeSearch(&(*(tmpKeys.begin())),&(*(pBalanced.begin())),0,tmpKeys.size(),0,pBalanced.size(),m_uiMaxDepth,m_uiMaxDepth,ROOT_ROTATION);
        SFC::seqSort::SFC_treeSort(&(*(tmpKeys.begin())),tmpKeys.size(),tmpNodes,tmpNodes,tmpNodes,m_uiMaxDepth,m_uiMaxDepth,root,ROOT_ROTATION,1,TS_SORT_ONLY);

        unsigned int keyCount=0;
        unsigned int nodeCount=0;

        while(keyCount<tmpKeys.size() && nodeCount<pBalanced.size())
        {
            if(tmpKeys[keyCount]<pBalanced[nodeCount])
                keyCount++;
            else if(tmpKeys[keyCount]>pBalanced[nodeCount])
                nodeCount++;
            else
                {
                    tmpKeyResults.push_back(keyCount);
                    nodeCount++;
                    keyCount++;

                }
        }


        double sfc_time=0.0;

        for(unsigned int iter=0;iter<NUM_ITERARTIONS;iter++)
        {
            tmpNodes.clear();
            tmpKeys.clear();
            tmpKeyResults.clear();
            genGauss(0.1,numKeys,dim,pts_keys);
            pts2Octants(tmpNodes,pts_keys,numKeys*dim,dim,m_uiMaxDepth);

            for(unsigned int i=0;i<tmpNodes.size();i++)
                tmpKeys.push_back(ot::TreeNode(tmpNodes[i]));


            t1=MPI_Wtime();
            //SFC::seqSearch::SFC_treeSearch(&(*(tmpKeys.begin())),&(*(pBalanced.begin())),0,tmpKeys.size(),0,pBalanced.size(),m_uiMaxDepth,m_uiMaxDepth,ROOT_ROTATION);
            SFC::seqSort::SFC_treeSort(&(*(tmpKeys.begin())),tmpKeys.size(),tmpNodes,tmpNodes,tmpNodes,m_uiMaxDepth,m_uiMaxDepth,root,ROOT_ROTATION,1,TS_SORT_ONLY);

            unsigned int keyCount=0;
            unsigned int nodeCount=0;

            while(keyCount<tmpKeys.size() && nodeCount<pBalanced.size())
            {
                if(tmpKeys[keyCount]<pBalanced[nodeCount])
                    keyCount++;
                else if(tmpKeys[keyCount]>pBalanced[nodeCount])
                    nodeCount++;
                else
                {
                    tmpKeyResults.push_back(keyCount);
                    nodeCount++;
                    keyCount++;

                }
            }

            t2=MPI_Wtime();
            sfc_time+=(t2-t1);
        }


        if(!rank) std::cout<<"sfcSearch(s) : "<<sfc_time<<"\t averge over "<<NUM_ITERARTIONS<<"(s): "<<(sfc_time/NUM_ITERARTIONS)<<std::endl;
        if(!rank) std::cout<<"number of keys: "<<tmpKeys.size()<<std::endl;
        //if(!rank)std::cout<<"last search found keys: "<<tmpKeyResults.size()<<std::endl;

        std::vector<unsigned int > foundKeys;
        double bsearch_time=0.0;

        // warmup run for bsearch
        foundKeys.clear();
        foundKeys.resize(numKeys,numKeys);
        tmpNodes.clear();
        tmpKeys.clear();
        genGauss(0.1,numKeys,dim,pts_keys);
        pts2Octants(tmpNodes,pts_keys,numKeys*dim,dim,m_uiMaxDepth);

        for(unsigned int i=0;i<tmpNodes.size();i++)
            tmpKeys.push_back(ot::Key(tmpNodes[i]));

        for(unsigned int i=0;i<tmpKeys.size();i++)
        {
            ot::TreeNode* found=(ot::TreeNode*)std::bsearch((ot::TreeNode*)&tmpKeys[i],&(*(pBalanced.begin())),pBalanced.size(),sizeof(ot::TreeNode),SFC::compare<ot::TreeNode>);
            if(found!=NULL)
                foundKeys[i]=((found-&(*(pBalanced.begin()))));
        }


        for(unsigned int iter=0;iter<NUM_ITERARTIONS;iter++)
        {

            foundKeys.clear();
            foundKeys.resize(numKeys,numKeys);
            tmpNodes.clear();
            tmpKeys.clear();
            genGauss(0.1,numKeys,dim,pts_keys);
            pts2Octants(tmpNodes,pts_keys,numKeys*dim,dim,m_uiMaxDepth);

            for(unsigned int i=0;i<tmpNodes.size();i++)
                tmpKeys.push_back(ot::Key(tmpNodes[i]));

            t1=MPI_Wtime();
            for(unsigned int i=0;i<tmpKeys.size();i++)
            {
                ot::TreeNode* found=(ot::TreeNode*)std::bsearch((ot::TreeNode*)&tmpKeys[i],&(*(pBalanced.begin())),pBalanced.size(),sizeof(ot::TreeNode),SFC::compare<ot::TreeNode>);
                if(found!=NULL)
                    foundKeys[i]=((found-&(*(pBalanced.begin()))));
            }
            t2=MPI_Wtime();
            bsearch_time+=(t2-t1);


        }



        if(!rank) std::cout<<"bsearch(s) : "<<bsearch_time<<"\t averge over "<<NUM_ITERARTIONS<<"(s): "<<(bsearch_time/NUM_ITERARTIONS)<<std::endl;
        if(!rank) std::cout<<"number of keys: "<<tmpKeys.size()<<std::endl;

        delete [] pts_keys;

    }

    pSorted.clear();
    pBalanced.clear();
    pConstructed.clear();


    MPI_Finalize();
    return 0;


}