//
// Created by milinda on 6/12/16.
//

/*
 * @author: Milinda Shayamal Fernando
 * School of Computing, University of Utah
 * @date: 6/12/2016
 *
 * This file contains code to test octree related utilities, like octree construction, balancing, implementations based on the treeSort aproach.
 *
 * */


#include "testOctUtils.h"


int readPtsFromFile(char* filename, std::vector<double>& pts) {
    FILE* infile;
    size_t res;
    infile = fopen(filename,"rb");
    unsigned int temp;
    res = fread(&temp,sizeof(unsigned int),1,infile);

    double* ptsTemp = NULL;

    if(temp) {
        ptsTemp = new double[3*temp];
        assert(ptsTemp);
        res = fread(ptsTemp, sizeof(double),3*temp,infile);
    }

    fclose(infile);
    pts.resize(3*temp);

    for (int i=0; i < (3*temp); i++) {
        pts[i] = ptsTemp[i];
    }//end for

    if(ptsTemp) {
        delete [] ptsTemp;
        ptsTemp = NULL;
    }

    // std::cout << __func__ << ": size " << temp << ", " << pts.size() << std::endl;
    return 1;
}//end function







int main(int argc, char** argv) {

    MPI_Init(&argc, &argv);

    int rank, npes;
    MPI_Comm GLOBAL_COMM=MPI_COMM_WORLD;
    MPI_Comm_rank(GLOBAL_COMM, &rank);
    MPI_Comm_size(GLOBAL_COMM, &npes);

    unsigned int options=0;


    if(argc<3)
    {   if(!rank)
            std::cout<<"Usage :"<<argv[0]<<" numPts "<<" dim "<<" maxDepth "<<" tol  "<<" distribution (0- Normal 1- Uniform 2- LogarithmicNormal) options( 1-Remove duplicates 2- constructOctree 6-balancedOctree)"<<std::endl;
    }

    DendroIntL numPts=atoll(argv[1]);
    unsigned int dim=atoi(argv[2]);
    unsigned int maxDepth=atoi(argv[3]);
    double tol=0.001;
    unsigned int distribution=0;
    if(argc>4)
        tol=atof(argv[4]);

    if(argc>5)
        distribution=atoi(argv[5]);

    if(argc>6)
        options=atoi(argv[6]);


     _InitializeHcurve(dim);
     if(!rank) std::cout << "Initialized Hcurves for dimention "<<dim << std::endl;

    DendroIntL localSz, totalSz;


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

     if (!rank) std::cout << "finished generating "<<dim<<" dimensional points" << std::endl;

    //Read pts from files
    /*if (!rank) {
        std::cout << RED " Reading  " << argv[1] << NRM << std::endl; // Point size
    }

    readPtsFromFile(FileName, pts);

    if (!rank) {
        std::cout << GRN " Finished reading  " << argv[1] << NRM << std::endl; // Point size
    }




    for (DendroIntL i = 0; i < totPts; i += 3) {
        if ((pts[i] > 0.0) &&
            (pts[i + 1] > 0.0)
            && (pts[i + 2] > 0.0) &&
            (((unsigned int) (pts[i] * ((double) (1u << maxDepth)))) < (1u << maxDepth)) &&
            (((unsigned int) (pts[i + 1] * ((double) (1u << maxDepth)))) < (1u << maxDepth)) &&
            (((unsigned int) (pts[i + 2] * ((double) (1u << maxDepth)))) < (1u << maxDepth))) {

            tmpNodes.push_back(ot::TreeNode((unsigned int) (pts[i] * (double) (1u << maxDepth)),
                                            (unsigned int) (pts[i + 1] * (double) (1u << maxDepth)),
                                            (unsigned int) (pts[i + 2] * (double) (1u << maxDepth)),
                                            maxDepth, dim, maxDepth));
        }
    }*/



     pts2Octants(tmpNodes,pts,totPts,dim,maxDepth); // generate octants from the points.

    //std::cout<<"inp size: "<<tmpNodes.size()<<std::endl;

    std::vector<ot::TreeNode> pNodesSorted;
    std::vector<ot::TreeNode> pNodesConstructed;
    std::vector<ot::TreeNode> pNodesBalanced;

    ot::TreeNode root(0,0,0,0,dim,maxDepth);


    localSz=tmpNodes.size();
    par::Mpi_Reduce(&localSz,&totalSz,1,MPI_SUM,ROOT_PROC,GLOBAL_COMM);

    if(!rank) std::cout<<YLW<<"================================================================="<<NRM<<std::endl;
    if(!rank) std::cout<<YLW<<" Total number of Points: "<<totalSz<<NRM<<std::endl;
    if(!rank) std::cout<<YLW<<"================================================================="<<NRM<<std::endl;


    if(npes==1)
    {

        if(!rank) std::cout <<RED<<"==================================================================================================="<<NRM<<std::endl;
        if(!rank) std::cout <<RED<<"           Executing the SEQUENTIAL TREESORT with options : "<<options<<NRM<<std::endl;
        if(!rank) std::cout<<RED<<"==================================================================================================="<<NRM<<std::endl;
        SFC::seqSort::SFC_treeSort(&(*(tmpNodes.begin())),tmpNodes.size(),pNodesSorted,pNodesConstructed,pNodesBalanced,maxDepth,maxDepth,root,0,1,options);
        assert(seq::test::isUniqueAndSorted(pNodesSorted));
        assert(seq::test::isUniqueAndSorted(pNodesConstructed));
        assert(seq::test::isUniqueAndSorted(pNodesBalanced));

        /*bool isUniqueSorted=seq::test::isUniqueAndSorted(pNodesSorted);

        if(isUniqueSorted)
        {
            if(!rank) std::cout<<"seq::test::isUniqueAndSorted ............. PASSED "<<std::endl;

        }else
        {
            if(!rank) std::cout<<"seq::test::isUniqueAndSorted ............. FAILED "<<std::endl;
        }




        bool isUniqueSorted= seq::test::isUniqueAndSorted(pNodesConstructed);
        bool isComplete=seq::test::isComplete(pNodesConstructed);

        if(isComplete & isUniqueSorted)
        {
            if(!rank) std::cout<<"seq::test::isCompleteOctree and seq::test::isUniqueAndSorted ............. PASSED "<<std::endl;
        }else
        {
            if(!rank) std::cout<<"seq::test::isCompleteOctree and seq::test::isUniqueAndSorted ............. FAILED "<<std::endl;
        }


        assert(isUniqueSorted);
        assert(isComplete);*/

    }else
    {

        if(!rank) std::cout <<RED<<"==================================================================================================="<<NRM<<std::endl;
        if(!rank) std::cout <<RED<<"           Executing the PARALLEL TREESORT with options : "<<options<<NRM<<std::endl;
        if(!rank) std::cout<<RED<<"==================================================================================================="<<NRM<<std::endl;

        SFC::parSort::SFC_treeSort(tmpNodes,pNodesSorted,pNodesConstructed,pNodesBalanced,tol,maxDepth,root,0,1,options,2,MPI_COMM_WORLD);
        assert(par::test::isUniqueAndSorted(pNodesSorted,MPI_COMM_WORLD));
        assert(par::test::isUniqueAndSorted(pNodesConstructed,MPI_COMM_WORLD));
        assert(par::test::isUniqueAndSorted(pNodesBalanced,MPI_COMM_WORLD));

        assert(par::test::containsAncestor(pNodesConstructed,MPI_COMM_WORLD));
        assert(par::test::containsAncestor(pNodesBalanced,MPI_COMM_WORLD));



    }

    localSz=pNodesSorted.size();
    par::Mpi_Reduce(&localSz,&totalSz,1,MPI_SUM,ROOT_PROC,GLOBAL_COMM);

    if(!rank) std::cout<<YLW<<"================================================================="<<NRM<<std::endl;
    if(!rank) std::cout<<YLW<<" Total number of sorted & unique points: "<<totalSz<<NRM<<std::endl;
    if(!rank) std::cout<<YLW<<"================================================================="<<NRM<<std::endl;


    localSz=pNodesConstructed.size();
    par::Mpi_Reduce(&localSz,&totalSz,1,MPI_SUM,ROOT_PROC,GLOBAL_COMM);

    if(!rank) std::cout<<YLW<<"================================================================="<<NRM<<std::endl;
    if(!rank) std::cout<<YLW<<" Total number of octants in Complete octree: "<<totalSz<<NRM<<std::endl;
    if(!rank) std::cout<<YLW<<"================================================================="<<NRM<<std::endl;


    localSz=pNodesBalanced.size();
    par::Mpi_Reduce(&localSz,&totalSz,1,MPI_SUM,ROOT_PROC,GLOBAL_COMM);

    if(!rank) std::cout<<YLW<<"================================================================="<<NRM<<std::endl;
    if(!rank) std::cout<<YLW<<" Total number of balanced octants : "<<totalSz<<NRM<<std::endl;
    if(!rank) std::cout<<YLW<<"================================================================="<<NRM<<std::endl;


    treeNodesTovtk(pNodesConstructed,rank,"completeOctree");
    treeNodesTovtk(pNodesBalanced,rank,"balOctree");
    MPI_Finalize();
    return 0;



}


