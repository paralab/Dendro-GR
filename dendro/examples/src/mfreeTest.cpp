//
// Created by milinda on 1/9/19.
//

#include "mfreeTest.h"

namespace mfree
{

    profiler_t tLev[31];


    void ptBucketRoofline(const Point* ptList, unsigned int numPts,unsigned int lev,unsigned int pMaxDepth)
    {

        if(numPts<=1 || lev<1) return;

        unsigned int counts[NUM_CHILDREN]={0};

        Point** ptChild=new Point*[NUM_CHILDREN];
        if(lev==pMaxDepth)
        {
            for(unsigned int l=0;l<pMaxDepth;l++)
                tLev[l].snapreset();
        }


        tLev[(pMaxDepth-lev)].start();


        #pragma omp parallel
        {
            unsigned int countsPrivate[NUM_CHILDREN]={0};

            #pragma omp for
            for(unsigned int i=0;i<numPts;i++)
                countsPrivate[i%NUM_CHILDREN]++;


            #pragma omp critical
            {
                for(unsigned int child=0;child<NUM_CHILDREN ;child++)
                    counts[child]+=countsPrivate[child];
            }

        }


        for(unsigned int child=0;child<NUM_CHILDREN;child++)
        {
            (counts[child]>0 ) ? ptChild[child]=new Point[counts[child]] : ptChild[child]=NULL;
        }

        for(unsigned int i=0;i<NUM_CHILDREN;i++)
            counts[i]=0;


        #pragma omp parallel
        {
            unsigned int cnum;

            #pragma omp for
            for(unsigned int i=0;i<numPts;i++)
            {
                cnum=(i%NUM_CHILDREN);
                ptChild[cnum][counts[cnum]]=ptList[i];
                #pragma omp critical
                {
                    counts[cnum]++;
                }

            }

        }

        tLev[(pMaxDepth-lev)].stop();

        for(unsigned int child=0;child<NUM_CHILDREN;child++)
        {
            if(counts[child] && lev>=1)
                ptBucketRoofline(ptChild[child],counts[child],lev-1,pMaxDepth);
        }



        if(lev==pMaxDepth)
        {
            for(unsigned int l=0;l<pMaxDepth;l++)
                std::cout<<"lev: "<<(l)<<" time (s) : "<<tLev[l].snap<<std::endl;
        }


        for(unsigned int child=0;child<NUM_CHILDREN;child++)
            delete [] ptChild[child];

        delete [] ptChild;





    }
}




int main(int argc, char** argv)
{


    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;

    int rank, npes;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    if(npes>1)
    {
        std::cout<<" "<<std::endl;
    }


    if (argc < 3) {
        if (!rank)
            std::cout << "Usage: " << argv[0]
                      << " maxDepth numpts eleOrder"
                      << std::endl;
        return 0;
    }

    m_uiMaxDepth = atoi(argv[1]);
    unsigned int numPts =atoi(argv[2]);
    unsigned int eOrder = atoi(argv[3]);


    double tBegin = 0, tEnd = 10, th = 0.01;

    if (!rank) {
        std::cout << YLW << "maxDepth: " << m_uiMaxDepth << NRM << std::endl;
        std::cout << YLW << "numPts: " << numPts << NRM << std::endl;
        std::cout << YLW << "eleOrder: " << eOrder << NRM << std::endl;

    }

    _InitializeHcurve(m_uiDim);

    const unsigned int totPts=numPts*m_uiDim;
    const double partition_tol=0.1;
    const unsigned int grainSz=100;

    std::vector<ot::TreeNode> tmpNodes;
    double * pts=new double[totPts];
    genGauss(0.15,numPts,m_uiDim,pts);  // Default case.
    pts2Octants(tmpNodes,pts,totPts,m_uiDim,m_uiMaxDepth); // generate octants from the points.

    ot::TreeNode root(m_uiDim,m_uiMaxDepth);
    std::vector<ot::TreeNode> tmpVec;
    double t1=MPI_Wtime();

    SFC::parSort::SFC_treeSort(tmpNodes,tmpVec,tmpVec,tmpVec,partition_tol,m_uiMaxDepth,root,ROOT_ROTATION,1,TS_REMOVE_DUPLICATES,DENDRO_DEFAULT_SF_K,comm);
    std::swap(tmpNodes,tmpVec);
    tmpVec.clear();

    SFC::parSort::SFC_treeSort(tmpNodes,tmpVec,tmpVec,tmpVec,partition_tol,m_uiMaxDepth,root,ROOT_ROTATION,1,TS_CONSTRUCT_OCTREE,DENDRO_DEFAULT_SF_K,comm);
    std::swap(tmpNodes,tmpVec);
    tmpVec.clear();

    double t2=MPI_Wtime();
    double t_stat=t2-t1;

    double t_stat_g[3];

    par::Mpi_Reduce(&t_stat,t_stat_g,1,MPI_MIN,0,comm);
    par::Mpi_Reduce(&t_stat,t_stat_g+1,1,MPI_SUM,0,comm);
    par::Mpi_Reduce(&t_stat,t_stat_g+2,1,MPI_MAX,0,comm);
    t_stat_g[1]=t_stat_g[1]/(double)npes;

    double localSz=tmpNodes.size();
    double globalSz;
    par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);

    if(!rank) std::cout<<GRN<<"remove duplicates + octree construction (s): "<<t_stat_g[2]<<NRM<<std::endl;
    if(!rank) std::cout<<GRN<<" # const. octants: "<<globalSz<<NRM<<std::endl;

    std::vector<ot::TreeNode> balOct;

    t1=MPI_Wtime();

    SFC::parSort::SFC_treeSort(tmpNodes,balOct,balOct,balOct,partition_tol,m_uiMaxDepth,root,ROOT_ROTATION,1,TS_BALANCE_OCTREE,DENDRO_DEFAULT_SF_K,comm);
    tmpNodes.clear();

    t2=MPI_Wtime();

    t_stat=t2-t1;
    par::Mpi_Reduce(&t_stat,t_stat_g,1,MPI_MIN,0,comm);
    par::Mpi_Reduce(&t_stat,t_stat_g+1,1,MPI_SUM,0,comm);
    par::Mpi_Reduce(&t_stat,t_stat_g+2,1,MPI_MAX,0,comm);
    t_stat_g[1]=t_stat_g[1]/(double)npes;

    if(!rank) std::cout<<GRN<<" 2:1 balancing max (s): "<<t_stat_g[2]<<NRM<<std::endl;
    localSz=balOct.size();


    par::Mpi_Allreduce(&localSz,&globalSz,1,MPI_SUM,comm);
    if(!rank) std::cout<<GRN<<" balanced # octants : "<<globalSz<<NRM<<std::endl;


    t1=MPI_Wtime();
    ot::Mesh * mesh=new ot::Mesh(balOct,1,eOrder,comm,true,ot::SM_TYPE::FEM_CG,grainSz,partition_tol);
    t2=MPI_Wtime();

    t_stat=t2-t1;
    par::Mpi_Reduce(&t_stat,t_stat_g,1,MPI_MIN,0,comm);
    par::Mpi_Reduce(&t_stat,t_stat_g+1,1,MPI_SUM,0,comm);
    par::Mpi_Reduce(&t_stat,t_stat_g+2,1,MPI_MAX,0,comm);
    t_stat_g[1]=t_stat_g[1]/(double)npes;

    localSz=mesh->getNumLocalMeshNodes();
    par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);
    if(!rank) std::cout<<GRN<<" # of CG nodes (vertices) : "<<globalSz<<NRM<<std::endl;
    if(!rank)
    {
        std::cout<< GRN<<"Mesh generation time (max): "<<t_stat_g[2]<<NRM<<std::endl;
        std::cout<<"\t"<<GRN<<" e2e (min,mean,max): "<<"( "<<t_e2e_g[0]<<"\t"<<t_e2e_g[1]<<"\t"<<t_e2e_g[2]<<" )"<<NRM<<std::endl;
        std::cout<<"\t"<<GRN<<" e2n (min,mean,max): "<<"( "<<t_e2n_g[0]<<"\t"<<t_e2n_g[1]<<"\t"<<t_e2n_g[2]<<" )"<<NRM<<std::endl;
        std::cout<<"\t"<<GRN<<" sm (min,mean,max): "<<"( "<<t_sm_g[0]<<"\t"<<t_sm_g[1]<<"\t"<<t_sm_g[2]<<" )"<<NRM<<std::endl;
        std::cout<<"\t"<<GRN<<" blk (min,mean,max): "<<"( "<<t_blk_g[0]<<"\t"<<t_blk_g[1]<<"\t"<<t_blk_g[2]<<" )"<<NRM<<std::endl;
    }


    const unsigned int dof=mesh->getDegOfFreedom();
    const unsigned int nodeLocalBegin=mesh->getNodeLocalBegin();
    const unsigned int nodeLocalEnd=mesh->getNodeLocalEnd();
    const unsigned int elementLocalBegin=mesh->getElementLocalBegin();
    const unsigned int elementLocalEnd=mesh->getElementLocalEnd();
    const unsigned int local_dof=nodeLocalEnd-nodeLocalBegin;

    ot::TreeNode rootNode(0,0,0,0,m_uiDim,m_uiMaxDepth);

    const unsigned int nPe=mesh->getNumNodesPerElement();
    const ot::TreeNode * allElements=&(*(mesh->getAllElements().begin()));
    const unsigned int * e2n=&(*(mesh->getE2NMapping().begin()));
    const unsigned int * e2n_dg=&(*(mesh->getE2NMapping_DG().begin()));
    const unsigned int * cg2dg=&(*(mesh->getCG2DGMap().begin()));
    const unsigned int eleOrder=mesh->getElementOrder();
    unsigned int nodeLookUp_CG,nodeLookUp_DG,ownerID,ii_x,jj_y,kk_z,len;
    unsigned int x,y,z;
    std::vector<unsigned int> cgNodeIndex;
    cgNodeIndex.resize(nodeLocalEnd-nodeLocalBegin);

    std::vector<Point> ptList;


    for(unsigned int node=nodeLocalBegin;node<nodeLocalEnd;node++)
    {
        nodeLookUp_DG=cg2dg[node];
        mesh->dg2eijk(nodeLookUp_DG,ownerID,ii_x,jj_y,kk_z);
        len=1u<<(m_uiMaxDepth-allElements[ownerID].getLevel());
        x=allElements[ownerID].getX()+ ii_x*(len/(eleOrder));
        y=allElements[ownerID].getY()+ jj_y*(len/(eleOrder));
        z=allElements[ownerID].getZ()+ kk_z*(len/(eleOrder));
        ptList.push_back(Point(x,y,z));
        cgNodeIndex[node-nodeLocalBegin]=node;

    }

    localSz=ptList.size();
    par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);



    if(!rank) std::cout<<"Warm up run"<<std::endl;
    mfree::ptBucketRoofline(&(*(ptList.begin())),ptList.size(),m_uiMaxDepth,m_uiMaxDepth);


    if(!rank) std::cout<<YLW<<"Total number of cg points :  "<<globalSz<<std::endl;
    profiler_t tBucket;
    tBucket.start();

    t1=MPI_Wtime();
    tBucket.start();
    mfree::ptBucketRoofline(&(*(ptList.begin())),ptList.size(),m_uiMaxDepth,m_uiMaxDepth);
    tBucket.stop();
    t2=MPI_Wtime();

    if(!rank) std::cout<<"SFC bucket time: "<<tBucket.snap<<" mpi: "<<(t2-t1)<<std::endl;

    double * vecIn;
    double * vecOut;
    vecIn=mesh->createVector<double>();
    vecOut=mesh->createVector<double>();

    #pragma omp parallel for
    for(unsigned int ele=elementLocalBegin;ele<elementLocalEnd;ele++)
    {
        for(unsigned int node=0;node<nPe;node++)
        {
            vecIn[e2n[ele*nPe+node]]=ele;
        }

    }

    profiler_t tIndirectMemAccess;
    tIndirectMemAccess.start();

    tIndirectMemAccess.start();
    #pragma omp parallel for
    for(unsigned int ele=elementLocalBegin;ele<elementLocalEnd;ele++)
    {
        for(unsigned int node=0;node<nPe;node++)
        {
            vecIn[e2n[ele*nPe+node]]=ele;
        }

    }

    tIndirectMemAccess.stop();


    if(!rank) std::cout<<"Indirect Memory Access time: "<<tIndirectMemAccess.snap<<std::endl;



    delete [] vecIn;
    delete [] vecOut;








    delete mesh;
    MPI_Finalize();
    return 0;

}