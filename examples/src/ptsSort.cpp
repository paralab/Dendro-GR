//
// Created by milinda on 2/11/16.
//

#include "iostream"
#include "TreeNode.h"
#include <vector>
#include <fstream>
#include "binUtils.h"
#include <algorithm>
#include <assert.h>
#include <testUtils.h>
#include "hcurvedata.h"
#include "testUtils.h"


int main(int  argc,char ** argv)
{

    if(argc<2)
    {
        std::cout<<"Usage :"<<argv[0]<<" ptsFilename"<<std::endl;
        return -1;
    }

    unsigned int numPts;

    _InitializeHcurve(3);

    double pts[3];
    char other[256];
    double xyzmax[3];

    std::ifstream myfile;
    myfile.open (argv[1]);

    if(!myfile.is_open())
    {
        std::cout<<"File open failed. "<<std::endl;
    }

    myfile>>numPts;

    myfile>>other;
    myfile>>other;

    myfile>>other;
    xyzmax[0]=atof(other);
    myfile>>other;
    xyzmax[1]=atof(other);
    myfile>>other;
    xyzmax[2]=atof(other);

    int maxLimit=std::max(xyzmax[0],std::max(xyzmax[1],xyzmax[2]));
    unsigned int depth=binOp::binLength(maxLimit+1);

    std::vector<ot::TreeNode> nodes;
    for(int i=0;i<numPts;i++)
    {
        myfile>>other;

        myfile>>other;
        pts[0]=atof(other);
        myfile>>other;
        pts[1]=atof(other);
        myfile>>other;
        pts[2]=atof(other);

        ot::TreeNode tmp(pts[0],pts[1],pts[2],depth,3,depth);
        nodes.push_back(tmp);
        //std::cout<<"Read Point: "<<tmp<<std::endl;
    }


    //std::cout<<"File Read Complete"<<std::endl;
    myfile.close();
    std::sort(nodes.begin(),nodes.end());
    assert(seq::test::isSorted(nodes));
    char ptsFileName[256];

#ifdef HILBERT_ORDERING
    sprintf(ptsFileName,"%s_%s",argv[1],"_H.xyz");
#else
    sprintf(ptsFileName,"%s_%s",argv[1],"_M.xyz");
#endif

    std::ofstream sorted_out;
    sorted_out.open(ptsFileName);
    sorted_out<<numPts<<std::endl;

    for(int i=0;i<numPts;i++)
        sorted_out<<"X\t"<<nodes[i].getX()<<"\t"<<nodes[i].getY()<<"\t"<<nodes[i].getZ()<<std::endl;

    sorted_out.close();








}


