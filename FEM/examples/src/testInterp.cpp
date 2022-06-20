//
// Created by milinda on 3/31/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief Creates a unit cube (element) and perform the interpolation by refining it hierarchically.
*/
//

#include "TreeNode.h"
#include "mpi.h"
#include "genPts_par.h"
#include "sfcSort.h"
#include "mesh.h"
#include "dendro.h"
#include "dendroIO.h"
#include "octUtils.h"
#include "functional"
int main (int argc, char** argv) {


    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;

    int rank, npes;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    if(argc<5)
    {
        if(!rank) std::cout<<"Usage: "<<argv[0]<<" maxDepth wavelet_tol partition_tol eleOrder cnum"<<std::endl;
        return 0;
    }

    m_uiMaxDepth=atoi(argv[1]);
    double wavelet_tol=atof(argv[2]);
    double partition_tol=atof(argv[3]);
    unsigned int eOrder=atoi(argv[4]);
    unsigned int cnum=atoi(argv[5]);

    if(!rank)
    {
        std::cout<<YLW<<"maxDepth: "<<m_uiMaxDepth<<NRM<<std::endl;
        std::cout<<YLW<<"wavelet_tol: "<<wavelet_tol<<NRM<<std::endl;
        std::cout<<YLW<<"partition_tol: "<<partition_tol<<NRM<<std::endl;
        std::cout<<YLW<<"eleOrder: "<<eOrder<<NRM<<std::endl;
        std::cout<<YLW<<"cnum: "<<cnum<<NRM<<std::endl;

    }

    _InitializeHcurve(m_uiDim);

    // function that we need to interpolate.
    std::function<double(double,double,double)> func =[](const double x,const double y,const double z){ return (sin(2*M_PI*x)*sin(2*M_PI*y)*sin(2*M_PI*z));};
    std::function<double(double)> func1D =[](const double x){ return (sin(2*M_PI*x));};
    //std::function<double(double,double,double)> func =[](const double &x,const double & y,const double &z){ return (x*x+y*y+z*z);};
    //std::function<double(double,double,double)> func =[](const double &x,const double & y,const double &z){ return (x+y+z);};

    std::vector<ot::TreeNode> tmpNodes;
    if(npes>1) if(!rank) std::cout<<"This needs to be run sequentially. "<<std::endl;

    if(npes==1)
    {
        //ot::TreeNode elem=ot::TreeNode(m_uiDim,m_uiMaxDepth);
        /*unsigned int dim=1;

        unsigned int nodesPerElem=(eOrder+1);
        double * pValue=new double[nodesPerElem];
        double * cValue=new double[nodesPerElem];
        double * cValue_ip=new double[nodesPerElem];

        double domain_min=0.15;
        double domain_max=0.3;
        double hp=(domain_max-domain_min)/eOrder;
        double mySz=(domain_max-domain_min);
        RefElement refEl(m_uiDim,eOrder);


        unsigned int depth=0;
        double l2=0;
        Point ptParent;
        Point ptChild;
        ptParent=Point(domain_min,0.0,0.0);
        do
        {

            //std::cout<<"parent_pt: "<<ptParent.x()<<" ,"<<ptParent.y()<<" ,"<<ptParent.z()<<std::endl;
              for(unsigned int i=0;i<(eOrder+1);i++)
                   pValue[i]=func1D(ptParent.x()+i*hp);

            std::cout<<" Parent: ";printArray_1D(pValue,nodesPerElem);

            refEl.I1D_Parent2Child(pValue,cValue_ip,cnum);
            hp=hp/2.0;
            mySz=mySz/2.0;

            if(cnum==0)
                ptChild=Point(ptParent.x(),ptParent.y(),ptParent.z());
            else if(cnum==1)
                ptChild=Point(ptParent.x()+mySz,ptParent.y(),ptParent.z());

            //std::cout<<"child_pt: "<<ptChild.x()<<" ,"<<ptChild.y()<<" ,"<<ptChild.z()<<std::endl;

            ptParent=Point(ptChild.x(),ptChild.y(),ptChild.z());

            for(unsigned int i=0;i<(eOrder+1);i++)
                  cValue[i]=func1D(ptChild.x()+i*hp);

            std::cout<<" child_ip: ";printArray_1D(cValue_ip,nodesPerElem);
            std::cout<<" child: ";printArray_1D(cValue,nodesPerElem);
            l2=normL2(cValue,cValue_ip,nodesPerElem)/normL2(cValue,nodesPerElem);
            depth++;
            std::cout<<"depth: "<<depth<<" relErr: "<<l2<<std::endl;


        }while(depth<m_uiMaxDepth && l2> wavelet_tol);*/




         unsigned int nodesPerElem=(eOrder+1)*(eOrder+1)*(eOrder+1);
        double * pValue=new double[nodesPerElem];
        double * cValue=new double[nodesPerElem];
        double * cValue_ip=new double[nodesPerElem];

        double domain_min=0.15;
        double domain_max=0.3;
        double hp=(domain_max-domain_min)/eOrder;
        double mySz=(domain_max-domain_min);
        std::cout<<"hp: "<<hp<<std::endl;

        RefElement refEl(m_uiDim,eOrder);
        unsigned int depth=0;
        double l2=0;
        Point ptParent;
        Point ptChild;
        ptParent=Point(domain_min,domain_min,domain_min);
        do
        {

            //std::cout<<"parent_pt: "<<ptParent.x()<<" ,"<<ptParent.y()<<" ,"<<ptParent.z()<<std::endl;

            for(unsigned int k=0;k<(eOrder+1);k++)
                for(unsigned int j=0;j<(eOrder+1);j++)
                    for(unsigned int i=0;i<(eOrder+1);i++)
                        pValue[k*(eOrder+1)*(eOrder+1)+j*(eOrder+1)+i]=func(ptParent.x()+i*hp,ptParent.y()+j*hp,ptParent.z()+k*hp);

           /* std::cout<<" Parent: ";printArray_1D(pValue,nodesPerElem);*/
            refEl.I3D_Parent2Child(pValue,cValue_ip,cnum);
            hp=hp/2.0;
            mySz=mySz/2.0;

            if(cnum==0)
                ptChild=Point(ptParent.x(),ptParent.y(),ptParent.z());
            else if(cnum==1)
                ptChild=Point(ptParent.x()+mySz,ptParent.y(),ptParent.z());
            else if(cnum==2)
                ptChild=Point(ptParent.x(),ptParent.y()+mySz,ptParent.z());
            else if(cnum==3)
                ptChild=Point(ptParent.x()+mySz,ptParent.y()+mySz,ptParent.z());
            else if(cnum==4)
                ptChild=Point(ptParent.x(),ptParent.y(),ptParent.z()+mySz);
            else if(cnum==5)
                ptChild=Point(ptParent.x()+mySz,ptParent.y(),ptParent.z()+mySz);
            else if(cnum==6)
                ptChild=Point(ptParent.x(),ptParent.y()+mySz,ptParent.z()+mySz);
            else if(cnum==7)
                ptChild=Point(ptParent.x()+mySz,ptParent.y()+mySz,ptParent.z()+mySz);

            //std::cout<<"child_pt: "<<ptChild.x()<<" ,"<<ptChild.y()<<" ,"<<ptChild.z()<<std::endl;

            ptParent=Point(ptChild.x(),ptChild.y(),ptChild.z());

            for(unsigned int k=0;k<(eOrder+1);k++)
                for(unsigned int j=0;j<(eOrder+1);j++)
                    for(unsigned int i=0;i<(eOrder+1);i++)
                        cValue[k*(eOrder+1)*(eOrder+1)+j*(eOrder+1)+i]=func(ptChild.x()+i*hp,ptChild.y()+j*hp,ptChild.z()+k*hp);

           /* std::cout<<" child_ip: ";printArray_1D(cValue_ip,nodesPerElem);
            std::cout<<" child: ";printArray_1D(cValue,nodesPerElem);*/

            l2=normL2(cValue,cValue_ip,nodesPerElem)/normL2(cValue,nodesPerElem);
            depth++;
            std::cout<<"depth: "<<depth<<" relErr: "<<l2<<std::endl;


        }while(depth<m_uiMaxDepth && l2> wavelet_tol);





    }

    MPI_Finalize();

}