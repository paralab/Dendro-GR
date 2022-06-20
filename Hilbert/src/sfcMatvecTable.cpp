//
// Created by milinda on 8/6/18.
//
/**
 * @author Milinda Fernando
 * School of Computing, University of Utah
 * @brief computes the sfc matvec bucketing table.
 * */

#include <iostream>
#include "dendro.h"
#include <math.h>
#include <vector>
#include <algorithm>




int main(int argc, char **argv)
{


    unsigned int rows;
    unsigned int r_y,r_z;

    if(m_uiDim==3)
    {

        // index :
        // 0 - min
        // 1 - mid
        // 2 - max


        rows=pow(3,m_uiDim);
        r_z=m_uiDim*m_uiDim;
        r_y=m_uiDim;



        std::vector<unsigned int > x_min;
        std::vector<unsigned int > x_mid;
        std::vector<unsigned int > x_max;

        std::vector<unsigned int > y_min;
        std::vector<unsigned int > y_mid;
        std::vector<unsigned int > y_max;

        std::vector<unsigned int > z_min;
        std::vector<unsigned int > z_mid;
        std::vector<unsigned int > z_max;


        x_min.push_back(0);
        x_min.push_back(2);
        x_min.push_back(4);
        x_min.push_back(6);

        x_max.push_back(1);
        x_max.push_back(3);
        x_max.push_back(5);
        x_max.push_back(7);

        y_min.push_back(0);
        y_min.push_back(1);
        y_min.push_back(4);
        y_min.push_back(5);

        y_max.push_back(2);
        y_max.push_back(3);
        y_max.push_back(6);
        y_max.push_back(7);

        z_min.push_back(0);
        z_min.push_back(1);
        z_min.push_back(2);
        z_min.push_back(3);

        z_max.push_back(4);
        z_max.push_back(5);
        z_max.push_back(6);
        z_max.push_back(7);


        for(unsigned int child=0;child<NUM_CHILDREN;child++)
        {
            x_mid.push_back(child);
            y_mid.push_back(child);
            z_mid.push_back(child);

        }

        unsigned int sfcTable[rows*NUM_CHILDREN];

        for(unsigned int r=0;r<rows;r++)
            for(unsigned int child=0;child<NUM_CHILDREN;child++)
                sfcTable[r*NUM_CHILDREN+child]=UCHAR_MAX;


        std::vector<unsigned int> * _x[]={&x_min,&x_mid,&x_max};
        std::vector<unsigned int> * _y[]={&y_min,&y_mid,&y_max};
        std::vector<unsigned int> * _z[]={&z_min,&z_mid,&z_max};


        std::vector<unsigned int> intersec_zy;
        std::vector<unsigned int> intersec_zyx;


        for(unsigned int k=0;k<m_uiDim;k++)
            for(unsigned int j=0;j<m_uiDim;j++)
                for(unsigned int i=0;i<m_uiDim;i++)
                {
                    intersec_zy.clear();
                    intersec_zyx.clear();

                    std::set_intersection(_z[k]->begin(),_z[k]->end(),_y[j]->begin(),_y[j]->end(),back_inserter(intersec_zy));
                    std::sort(intersec_zy.begin(),intersec_zy.end());

                    std::set_intersection(intersec_zy.begin(),intersec_zy.end(),_x[i]->begin(),_x[i]->end(),back_inserter(intersec_zyx));
                    std::sort(intersec_zyx.begin(),intersec_zyx.end());

                    for(unsigned int ele=0;ele<intersec_zyx.size();ele++)
                        sfcTable[(k*r_z+j*r_y+i)*NUM_CHILDREN+ele]=intersec_zyx[ele];

                }



        std::cout<<"/**@brief generated sfc matvec bucking table, do not change !!! */ \n\n"<<std::endl;
        std::cout<<"//m_uiDim  "<<m_uiDim<<" \n\n"<<std::endl;
        std::cout<<"#define SFC_MVEC_TABLE_DEFAULT UCHAR_MAX"<<std::endl;
        std::cout<<"#define SFC_MVEC_TABLE_Rz 9"<<std::endl;
        std::cout<<"#define SFC_MVEC_TABLE_Ry 3"<<std::endl;

        std::cout<<"static const unsigned char SFC_MATVEC_BUCKET_TABLE[][]={ "<<std::endl;
        for(unsigned int r=0;r<rows;r++)
        {
            std::cout<<" {";
            for(unsigned int child=0;child<NUM_CHILDREN;child++)
                std::cout<<sfcTable[r*NUM_CHILDREN+child]<<",";

             std::cout<<"\b},"<<std::endl;
        }

        std::cout<<"};"<<std::endl;




    }








}

