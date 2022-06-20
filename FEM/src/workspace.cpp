//
// Created by milinda on 12/29/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief workspace variables allocations
*/
//

#include "workspace.h"

namespace fem
{
    namespace domain
    {
        Point grid_min;
        Point grid_max;

        Point domain_min;
        Point domain_max;

        double Rg_x;
        double Rg_y;
        double Rg_z;

        double Rd_x;
        double Rd_y;
        double Rd_z;

    }// end of namespace domain.

} // end of namespace fem



namespace fem
{
    namespace operators
    {

        namespace poisson
        {
            double *Qx;
            double *Qy;
            double *Qz;

            double *imV1;
            double *imV2;
            //double *imV3;

            double * pts_x;
            double * pts_y;
            double * pts_z;


            double ** jdXbydR;
            double ** jdRbydX;

            double *  jfactor;


            void allocateWorkSpace(const unsigned int n)
            {
                Qx=new double[n];
                Qy=new double[n];
                Qz=new double[n];

                imV1=new double[n];
                imV2=new double[n];
                //imV3=new double[n];

                jdXbydR=new double*[JACOBIAN_3D_SIZE];
                jdRbydX=new double*[JACOBIAN_3D_SIZE];

                for(unsigned int i=0;i<JACOBIAN_3D_SIZE;i++)
                {
                    jdXbydR[i]=new double [n];
                    jdRbydX[i]=new double [n];

                }

                jfactor=new double[n];

                pts_x=new double[n];
                pts_y=new double[n];
                pts_z=new double[n];



            }


            void deallocateWorkSpace()
            {
                delete [] Qx;
                delete [] Qy;
                delete [] Qz;

                delete [] imV1;
                delete [] imV2;
                //delete [] imV3;


                for(unsigned int i=0;i<JACOBIAN_3D_SIZE;i++)
                {
                    delete [] jdXbydR[i];
                    delete [] jdRbydX[i];

                }

                delete [] jdXbydR;
                delete [] jdRbydX;
                delete [] jfactor;

                delete [] pts_x;
                delete [] pts_y;
                delete [] pts_z;

            }

        }


    }// end of namespace operators

}// end of namespace fem


