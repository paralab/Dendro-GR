//
// Created by milinda on 1/23/17.
//

/**
 * @author: Milinda Fernanado
 * School of Computing, University of Utah.
 * @brief This an example of using the interpolation on a octree using tensor contraction.
 *
 * */



#include "TreeNode.h"
#include "sfcSort.h"
#include "mpi.h"
#include "hcurvedata.h"
#include "refel.h"
#include "genPts_par.h"

/*
#define DEBUG_INTERP_OCTREE_0
#define DEBUG_INTERP_OCTREE_1
#define DEBUG_INTERP_OCTREE_2
#define DEBUG_INTERP_OCTREE_3
#define DEBUG_INTERP_OCTREE_4
#define DEBUG_INTERP_OCTREE_5
#define DEBUG_INTERP_OCTREE_6
#define DEBUG_INTERP_OCTREE_7
*/


double evalOct(double x,double y, double z)
{
    //double val=x+y+z;
    double val=sin(x)+sin(y)+sin(z);
    //double val=x*x+y*y+z*z;
    //double val=sin(sqrt(x*x+y*y+z*z));//+y+z);
    return val;
}

double l2Norm(double * x, double* y, int n)
{
    double norm=0;
    for(unsigned int k=0;k<n;k++)
        norm+=((x[k]-y[k])*(x[k]-y[k]));

    //std::cout<<"norm: "<<norm<<std::endl;
    return sqrt(norm);
}


int main(int argc,char** argv)
{

    MPI_Init(&argc, &argv);
    int rank, npes;
    MPI_Comm GLOBAL_COMM=MPI_COMM_WORLD;
    MPI_Comm_rank(GLOBAL_COMM, &rank);
    MPI_Comm_size(GLOBAL_COMM, &npes);

    unsigned int options=0;


    if(argc<3)
    {   if(!rank)
            std::cout<<"Usage :"<<argv[0]<<" numPts "<<" dim "<<" maxDepth "<<" tol  "<<" distribution (0- Normal 1- Uniform 2- LogarithmicNormal) options( 1-Remove duplicates 2- constructOctree 6-balancedOctree)  splitterFixK  interpOrder"<<std::endl;
    }

    DendroIntL numPts=atoll(argv[1]);
    unsigned int dim=atoi(argv[2]);
    unsigned int maxDepth=atoi(argv[3]);
    m_uiMaxDepth=maxDepth;
    double tol=0.001;
    unsigned int distribution=0;
    unsigned int sf_k=2;
    unsigned int interpOrder=4;
    if(argc>4)
        tol=atof(argv[4]);

    if(argc>5)
        distribution=atoi(argv[5]);

    if(argc>6)
        options=atoi(argv[6]);

    if(argc>7)
        sf_k=atoi(argv[7]);

    if(argc>8)
        interpOrder=atoi(argv[8]);

    if(!rank) std::cout<<"sf parameter: "<<sf_k<<std::endl;

    _InitializeHcurve(dim);
    if(!rank) std::cout << "Initialized H-Curves for dimension "<<dim << std::endl;

    DendroIntL totPts=numPts*dim;
    double * pts=new double[totPts];

    std::vector<ot::TreeNode> tmpNodes;

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
    ot::TreeNode root(0,0,0,0,dim,maxDepth);


    if(npes>1) {
        SFC::parSort::SFC_treeSort(tmpNodes, pNodesSorted, pNodesConstructed, pNodesBalanced, tol, maxDepth, root,
                                   ROOT_ROTATION, 1, options, sf_k, GLOBAL_COMM);
    }else
    {
        SFC::seqSort::SFC_treeSort(&(*(tmpNodes.begin())),tmpNodes.size(),pNodesSorted,pNodesConstructed,pNodesBalanced,maxDepth,maxDepth,root,ROOT_ROTATION,1,options);
    }

    std::vector<ot::TreeNode> activeOctree;
    if(options==TS_REMOVE_DUPLICATES)
        std::swap(activeOctree,pNodesSorted);
    else if(options == TS_CONSTRUCT_OCTREE)
        std::swap(activeOctree,pNodesConstructed);
    else if(options == TS_BALANCE_OCTREE)
        std::swap(activeOctree,pNodesBalanced);

    MPI_Barrier(GLOBAL_COMM);
    if(!rank) std::cout<<"Octree Generated: "<<activeOctree.size()<<std::endl;

    RefElement refEl=RefElement(1,interpOrder);
    unsigned int Np=(interpOrder+1)*(interpOrder+1)*(interpOrder+1);
    double x,y, z;
    /*
     *
     * 1st Dim: x coords then y and then z
     *
     * */

    if(activeOctree.size()!=0)
    {
        /*activeOctree.clear();
        activeOctree.push_back(ot::TreeNode(0,0,0,0,m_uiDim,maxDepth));*/
        double* f_val=new double[activeOctree.size()*Np];
        for(unsigned int w=0;w<activeOctree.size();w++)
        {
            for(unsigned int k=0;k<(interpOrder+1);k++)
            {
                z=activeOctree[w].minZ()+ k*(activeOctree[w].maxZ()-activeOctree[w].minZ())/(double)interpOrder;

                for(unsigned int j=0;j<(interpOrder+1);j++)
                {
                    y=activeOctree[w].minY()+ j*(activeOctree[w].maxY()-activeOctree[w].minY())/(double)interpOrder;
                    for(unsigned int i=0;i<(interpOrder+1);i++)
                    {
                        x=activeOctree[w].minX()+ i*(activeOctree[w].maxX()-activeOctree[w].minX())/(double)interpOrder;
                        f_val[w*Np+k*(interpOrder+1)*(interpOrder+1)+j*(interpOrder+1)+i]=evalOct(x,y,z);
                    }

                }
            }

        }


        if(!rank) std::cout<<"Function Initialization Ended "<<std::endl;

        ot::TreeNode tmp;
        unsigned int xi,yi,zi,li,sz;
        double* f_calc=new double [Np];
        double* f_intp=new double [Np];

        for(unsigned int w=0;w<activeOctree.size();w++)
        {
            xi=activeOctree[w].getX();
            yi=activeOctree[w].getY();
            zi=activeOctree[w].getZ();
            li=activeOctree[w].getLevel();

            if(li<(m_uiMaxDepth-1)) {
                sz=1u<<(m_uiMaxDepth-li-1);
                for (unsigned int mChildNum = 0; mChildNum < NUM_CHILDREN; mChildNum++) {

                    switch (mChildNum) {
                        case 0:
                            tmp = ot::TreeNode(xi, yi, zi, li + 1, m_uiDim, m_uiMaxDepth);
                            refEl.I3D_Parent2Child(f_val+w*Np,f_intp,0);
                            for(unsigned int k=0;k<(interpOrder+1);k++) {
                                z = tmp.minZ() +
                                    k * (tmp.maxZ() - tmp.minZ()) / (double) interpOrder;
                                for (unsigned int j = 0; j < (interpOrder + 1); j++) {
                                    y = tmp.minY() +
                                        j * (tmp.maxY() - tmp.minY()) / (double) interpOrder;
                                    for (unsigned int i = 0; i < (interpOrder + 1); i++) {

                                        x = tmp.minX() +
                                            i * (tmp.maxX() - tmp.minX()) / (double) interpOrder;
                                        f_calc[k*(interpOrder+1)*(interpOrder+1)+j*(interpOrder+1)+i]=evalOct(x,y,z);

                                    }
                                }
                            }
#ifdef DEBUG_INTERP_OCTREE_0
                            if(!rank) std::cout<<" f_val for "<<w<<" :";printArray_1D(f_val+w*Np,Np);
                            if(!rank) std::cout<<" f_cal for "<<w<<" :";printArray_1D(f_calc,Np);
                            if(!rank) std::cout<<" f_inp for "<<w<<" :";printArray_1D(f_intp,Np);
#endif
                            break;

                        case 1:
                            tmp = ot::TreeNode(xi+sz, yi, zi, li + 1, m_uiDim, m_uiMaxDepth);
                            refEl.I3D_Parent2Child(f_val+w*Np,f_intp,1);
                            for(unsigned int k=0;k<(interpOrder+1);k++) {
                                z = tmp.minZ() +
                                    k * (tmp.maxZ() - tmp.minZ()) / (double) interpOrder;
                                for (unsigned int j = 0; j < (interpOrder + 1); j++) {
                                    y = tmp.minY() +
                                        j * (tmp.maxY() - tmp.minY()) / (double) interpOrder;
                                    for (unsigned int i = 0; i < (interpOrder + 1); i++) {

                                        x = tmp.minX() +
                                            i * (tmp.maxX() - tmp.minX()) / (double) interpOrder;
                                        f_calc[k*(interpOrder+1)*(interpOrder+1)+j*(interpOrder+1)+i]=evalOct(x,y,z);

                                    }
                                }
                            }
#ifdef DEBUG_INTERP_OCTREE_1
                            if(!rank) std::cout<<" f_val for "<<w<<" :";printArray_1D(f_val+w*Np,Np);
                            if(!rank) std::cout<<" f_cal for "<<w<<" :";printArray_1D(f_calc,Np);
                            if(!rank) std::cout<<" f_inp for "<<w<<" :";printArray_1D(f_intp,Np);
#endif
                            break;
                        case 2:
                            tmp = ot::TreeNode(xi, yi+sz, zi, li + 1, m_uiDim, m_uiMaxDepth);
                            refEl.I3D_Parent2Child(f_val+w*Np,f_intp,2);
                            for(unsigned int k=0;k<(interpOrder+1);k++) {
                                z = tmp.minZ() +
                                    k * (tmp.maxZ() - tmp.minZ()) / (double) interpOrder;
                                for (unsigned int j = 0; j < (interpOrder + 1); j++) {
                                    y = tmp.minY() +
                                        j * (tmp.maxY() - tmp.minY()) / (double) interpOrder;
                                    for (unsigned int i = 0; i < (interpOrder + 1); i++) {

                                        x = tmp.minX() +
                                            i * (tmp.maxX() - tmp.minX()) / (double) interpOrder;
                                        f_calc[k*(interpOrder+1)*(interpOrder+1)+j*(interpOrder+1)+i]=evalOct(x,y,z);

                                    }
                                }
                            }
#ifdef DEBUG_INTERP_OCTREE_2
                            if(!rank) std::cout<<" f_val for "<<w<<" :";printArray_1D(f_val+w*Np,Np);
                            if(!rank) std::cout<<" f_cal for "<<w<<" :";printArray_1D(f_calc,Np);
                            if(!rank) std::cout<<" f_inp for "<<w<<" :";printArray_1D(f_intp,Np);
#endif
                            break;

                        case 3:
                            tmp = ot::TreeNode(xi+sz, yi+sz, zi, li + 1, m_uiDim, m_uiMaxDepth);
                            refEl.I3D_Parent2Child(f_val+w*Np,f_intp,3);
                            for(unsigned int k=0;k<(interpOrder+1);k++) {
                                z = tmp.minZ() +
                                    k * (tmp.maxZ() - tmp.minZ()) / (double) interpOrder;
                                for (unsigned int j = 0; j < (interpOrder + 1); j++) {
                                    y = tmp.minY() +
                                        j * (tmp.maxY() - tmp.minY()) / (double) interpOrder;
                                    for (unsigned int i = 0; i < (interpOrder + 1); i++) {

                                        x = tmp.minX() +
                                            i * (tmp.maxX() - tmp.minX()) / (double) interpOrder;
                                        f_calc[k*(interpOrder+1)*(interpOrder+1)+j*(interpOrder+1)+i]=evalOct(x,y,z);

                                    }
                                }
                            }
#ifdef DEBUG_INTERP_OCTREE_3
                            if(!rank) std::cout<<" f_val for "<<w<<" :";printArray_1D(f_val+w*Np,Np);
                            if(!rank) std::cout<<" f_cal for "<<w<<" :";printArray_1D(f_calc,Np);
                            if(!rank) std::cout<<" f_inp for "<<w<<" :";printArray_1D(f_intp,Np);
#endif
                            break;
                        case 4:
                            tmp = ot::TreeNode(xi, yi, zi+sz, li + 1, m_uiDim, m_uiMaxDepth);
                            refEl.I3D_Parent2Child(f_val+w*Np,f_intp,4);
                            for(unsigned int k=0;k<(interpOrder+1);k++) {
                                z = tmp.minZ() +
                                    k * (tmp.maxZ() - tmp.minZ()) / (double) interpOrder;
                                for (unsigned int j = 0; j < (interpOrder + 1); j++) {
                                    y = tmp.minY() +
                                        j * (tmp.maxY() - tmp.minY()) / (double) interpOrder;
                                    for (unsigned int i = 0; i < (interpOrder + 1); i++) {

                                        x = tmp.minX() +
                                            i * (tmp.maxX() - tmp.minX()) / (double) interpOrder;
                                        f_calc[k*(interpOrder+1)*(interpOrder+1)+j*(interpOrder+1)+i]=evalOct(x,y,z);

                                    }
                                }
                            }
#ifdef DEBUG_INTERP_OCTREE_4
                            if(!rank) std::cout<<" f_val for "<<w<<" :";printArray_1D(f_val+w*Np,Np);
                            if(!rank) std::cout<<" f_cal for "<<w<<" :";printArray_1D(f_calc,Np);
                            if(!rank) std::cout<<" f_inp for "<<w<<" :";printArray_1D(f_intp,Np);
#endif
                            break;
                        case 5:
                            tmp = ot::TreeNode(xi+sz, yi, zi+sz, li + 1, m_uiDim, m_uiMaxDepth);
                            refEl.I3D_Parent2Child(f_val+w*Np,f_intp,5);
                            for(unsigned int k=0;k<(interpOrder+1);k++) {
                                z = tmp.minZ() +
                                    k * (tmp.maxZ() - tmp.minZ()) / (double) interpOrder;
                                for (unsigned int j = 0; j < (interpOrder + 1); j++) {
                                    y = tmp.minY() +
                                        j * (tmp.maxY() - tmp.minY()) / (double) interpOrder;
                                    for (unsigned int i = 0; i < (interpOrder + 1); i++) {

                                        x = tmp.minX() +
                                            i * (tmp.maxX() - tmp.minX()) / (double) interpOrder;
                                        f_calc[k*(interpOrder+1)*(interpOrder+1)+j*(interpOrder+1)+i]=evalOct(x,y,z);

                                    }
                                }
                            }
#ifdef DEBUG_INTERP_OCTREE_5
                            if(!rank) std::cout<<" f_val for "<<w<<" :";printArray_1D(f_val+w*Np,Np);
                            if(!rank) std::cout<<" f_cal for "<<w<<" :";printArray_1D(f_calc,Np);
                            if(!rank) std::cout<<" f_inp for "<<w<<" :";printArray_1D(f_intp,Np);
#endif
                            break;
                        case 6:
                            tmp = ot::TreeNode(xi, yi+sz, zi+sz, li + 1, m_uiDim, m_uiMaxDepth);
                            refEl.I3D_Parent2Child(f_val+w*Np,f_intp,6);
                            for(unsigned int k=0;k<(interpOrder+1);k++) {
                                z = tmp.minZ() +
                                    k * (tmp.maxZ() - tmp.minZ()) / (double) interpOrder;
                                for (unsigned int j = 0; j < (interpOrder + 1); j++) {
                                    y = tmp.minY() +
                                        j * (tmp.maxY() - tmp.minY()) / (double) interpOrder;
                                    for (unsigned int i = 0; i < (interpOrder + 1); i++) {

                                        x = tmp.minX() +
                                            i * (tmp.maxX() - tmp.minX()) / (double) interpOrder;
                                        f_calc[k*(interpOrder+1)*(interpOrder+1)+j*(interpOrder+1)+i]=evalOct(x,y,z);

                                    }
                                }
                            }
#ifdef DEBUG_INTERP_OCTREE_6
                            if(!rank) std::cout<<" f_val for "<<w<<" :";printArray_1D(f_val+w*Np,Np);
                            if(!rank) std::cout<<" f_cal for "<<w<<" :";printArray_1D(f_calc,Np);
                            if(!rank) std::cout<<" f_inp for "<<w<<" :";printArray_1D(f_intp,Np);
#endif
                            break;
                        case 7:
                            tmp = ot::TreeNode(xi+sz, yi+sz, zi+sz, li + 1, m_uiDim, m_uiMaxDepth);
                            refEl.I3D_Parent2Child(f_val+w*Np,f_intp,7);
                            for(unsigned int k=0;k<(interpOrder+1);k++) {
                                z = tmp.minZ() +
                                    k * (tmp.maxZ() - tmp.minZ()) / (double) interpOrder;
                                for (unsigned int j = 0; j < (interpOrder + 1); j++) {
                                    y = tmp.minY() +
                                        j * (tmp.maxY() - tmp.minY()) / (double) interpOrder;
                                    for (unsigned int i = 0; i < (interpOrder + 1); i++) {

                                        x = tmp.minX() +
                                            i * (tmp.maxX() - tmp.minX()) / (double) interpOrder;
                                        f_calc[k*(interpOrder+1)*(interpOrder+1)+j*(interpOrder+1)+i]=evalOct(x,y,z);

                                    }
                                }
                            }
#ifdef DEBUG_INTERP_OCTREE_7
                            if(!rank) std::cout<<" f_val for "<<w<<" :";printArray_1D(f_val+w*Np,Np);
                            if(!rank) std::cout<<" f_cal for "<<w<<" :";printArray_1D(f_calc,Np);
                            if(!rank) std::cout<<" f_inp for "<<w<<" :";printArray_1D(f_intp,Np);
#endif
                            break;

                    };




                    if(l2Norm(f_calc,f_intp,Np)>1) std::cout<<"Rank "<<rank<<" Octant: "<<w<<" to  Child Num : "<<mChildNum<<" normed difference: "<<l2Norm(f_calc,f_intp,Np)<<std::endl;



                }
            }
        }




        delete [] f_val;
        delete [] f_calc;
        delete [] f_intp;

    }








   MPI_Finalize();
   return 0;





}
