//
// Created by milinda on 1/20/17.
//

/**
 *
 * @author Milinda Fernando
 * School of Computing, University of Utah
 * @brief Contains example of generation refelement instance.
 * This tries to interpolate a sin(\theta)+cos(\theta) using legendre polynomials.
 *
 * */

#include "assert.h"
#include <iostream>
#include <math.h>
#include <cmath>
#include "cstring"

#include "refel.h"


double l2Norm(double * x, double* y, int n)
{
    double norm=0;
    for(unsigned int k=0;k<n;k++)
        norm+=((x[k]-y[k])*(x[k]-y[k]));

    //std::cout<<"norm: "<<norm<<std::endl;
    return sqrt(norm);
}


int main(int argc, char ** argv)
{
    if(argc<5)
    {
        std::cout<<"Usage: "<<argv[0]<<" dim "<<"order "<<"x_min "<<"x_max "<<" freq[enter k for k*pi]"<<std::endl;
        return 0;
    }

    int dim, max_order;
    dim=atoi(argv[1]);
    max_order=atoi(argv[2]);

    double x_min=atof(argv[3]);
    double x_max=atof(argv[4]);
    int freq_k=atoi(argv[5]);

    //RefElement rr=RefElement(1,10);

    double ** r =new double *[max_order];
    double ** f =new double *[max_order];

    double ** r0 =new double *[max_order];
    double ** f0 =new double *[max_order];
    double ** fI0 =new double *[max_order];


    double ** r1 =new double *[max_order];
    double ** f1 =new double *[max_order];
    double ** fI1 =new double *[max_order];




    for(unsigned int k=0;k<max_order;k++)
    {
        RefElement refEl=RefElement(dim,k+1);

        const double *I0=refEl.getIMChild0();
        const double *I1=refEl.getIMChild1();

        std::cout<<"child 0 : order: "<<(k+1)<<" I0: "<<std::endl;printArray_2D(I0,(k+2),(k+2));
        std::cout<<"child 1 : order: "<<(k+1)<<" I1: "<<std::endl;printArray_2D(I1,(k+2),(k+2));


        r[k]=new double[k+2];
        f[k]=new double[k+2];

        r0[k]=new double[k+2];
        f0[k]=new double[k+2];
        fI0[k]=new double[k+2];

        r1[k]=new double[k+2];
        f1[k]=new double[k+2];
        fI1[k]=new double[k+2];

        for(unsigned int w=0;w<(k+2);w++) {
            r[k][w] = x_min + w * (x_max - x_min) / (k + 1);
            f[k][w] = std::sin(k*M_PI*r[k][w])+std::cos(k*M_PI*r[k][w]);


            r0[k][w] = x_min + w * (0.5*(x_max+x_min) - x_min) / (k + 1);
            f0[k][w] = std::sin(k*M_PI*r0[k][w])+std::cos(k*M_PI*r0[k][w]);

            r1[k][w] = 0.5*(x_max+x_min) + w * (x_max - 0.5*(x_max+x_min)) / (k + 1);
            f1[k][w] = std::sin(k*M_PI*r1[k][w])+std::cos(k*M_PI*r1[k][w]);

            fI0[k][w]=0;
            fI1[k][w]=0;

        }




        for(unsigned int i=0;i<(k+2);i++)
        {
            for(unsigned int j=0;j<(k+2);j++)
            {
                //fI0[k][i]+=I0[i*(k+2)+j]*f0[k][j];
                //fI1[k][i]+=I1[i*(k+2)+j]*f1[k][j];

                fI0[k][i]+=f0[k][i]*I0[j*(k+2)+i];
                fI1[k][i]+=f1[k][i]*I1[j*(k+2)+i];
            }
        }

        std::cout<<"order: "<<(k+1)<<" r: "; printArray_1D(r[k],(k+2));
        std::cout<<"order: "<<(k+1)<<" f: "; printArray_1D(f[k],(k+2));

        std::cout<<"order: "<<(k+1)<<" r0: "; printArray_1D(r0[k],(k+2));
        std::cout<<"order: "<<(k+1)<<" f0: "; printArray_1D(f0[k],(k+2));
        std::cout<<"order: "<<(k+1)<<" fI0: "; printArray_1D(fI0[k],(k+2));

        std::cout<<"order: "<<(k+1)<<" r1: "; printArray_1D(r1[k],(k+2));
        std::cout<<"order: "<<(k+1)<<" f1: "; printArray_1D(f1[k],(k+2));
        std::cout<<"order: "<<(k+1)<<" fI1: "; printArray_1D(fI1[k],(k+2));

        std::cout<<"child 0 : order: "<<(k+1)<<" l2 norm(actual- interpolated): "<<l2Norm(f0[k],fI0[k],(k+2))<<std::endl;
        std::cout<<"child 1 : order: "<<(k+1)<<" l2 norm(actual- interpolated): "<<l2Norm(f1[k],fI1[k],(k+2))<<std::endl;




    }



    for(unsigned int i=0;i<max_order;i++) {
        delete[] r[i];
        delete[] f[i];

        delete[] r0[i];
        delete[] f0[i];

        delete[] r1[i];
        delete[] f1[i];

        delete[] fI0[i];
        delete[] fI1[i];
    }

    //delete [] refEl;




    return 0;


}