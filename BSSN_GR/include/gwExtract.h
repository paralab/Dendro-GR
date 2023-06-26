/**
 *@author: Milinda Fernando.
 *School of Computing University of Utah
 *@brief: Contains the COG/SymPyGR (symbolic) code generation to perform gravitational wave (GW)
 *extraction.
 *
 * The extraction method is based on "Extraction of Gravitational Waves in Numerical Relativity"
 * https://arxiv.org/abs/1606.02532 : Chapter 2 , 3, 4.
 *
 * For Spin Weighted Spherical Harmonics (SWSH) we use the SphericalFunction module based on
 * https://arxiv.org/pdf/1604.08140.pdf
 *
 *
 * */


#ifndef DENDRO_5_0_GWEXTRACT_H
#define DENDRO_5_0_GWEXTRACT_H

#include <iostream>
#include "grDef.h"
#include "mesh.h"
#include "parameters.h"
#include "daUtils.h"
#include <iomanip>
#include <complex>
#include <cmath>
#include "dendro.h"
#include "point.h"

extern "C" void dgesv_( int* n, int* nrhs, double* a, int* lda, int* ipiv,double* b, int* ldb, int* info );

namespace GW
{

    // Definitions related to SWHS and Lebedev quad. to perform GW extraction. 
    #include "nrswsh.h"

    /**
     * @brief extract the GW waves using poly fit.
     * @param [in] mesh: input mesh
     * @param [in] cVar: constraint variables.
     * @param [in] timestep : timestep value.
     * @param [in] time: simulation time.
     *
     * */

    template<typename T>
    void extractFarFieldPsi4(const ot::Mesh* mesh, const T** cVar,unsigned int timestep,double time)
    {
        // get the active comm and information about it
        unsigned int rankActive=mesh->getMPIRank();
        unsigned int npesActive=mesh->getMPICommSize();
        MPI_Comm commActive=mesh->getMPICommunicator();

        unsigned int totalModes=0;
        for(unsigned int l=0;l<BSSN_GW_NUM_LMODES;l++)
            totalModes+=2*BSSN_GW_L_MODES[l]+1;

        const unsigned int TOTAL_MODES=totalModes;

        DendroComplex * swsh_coeff = new DendroComplex[BSSN_GW_NUM_RADAII*TOTAL_MODES];
        DendroComplex * swsh_coeff_g = new DendroComplex[BSSN_GW_NUM_RADAII*TOTAL_MODES];

        std::vector<unsigned int> lmCounts;
        std::vector<unsigned int> lmOffset;

        lmCounts.resize(BSSN_GW_NUM_LMODES);
        lmOffset.resize(BSSN_GW_NUM_LMODES);

        for(unsigned int l=0;l<BSSN_GW_NUM_LMODES;l++)
            lmCounts[l]=2*BSSN_GW_L_MODES[l]+1;

        lmOffset[0]=0;
        omp_par::scan(&(*(lmCounts.begin())),&(*(lmOffset.begin())),BSSN_GW_NUM_LMODES);


        std::vector<double> psi4L2R;
        std::vector<double> psi4L2I;

        std::vector<double> psi4L2R_g;
        std::vector<double> psi4L2I_g;

        psi4L2R.resize(BSSN_GW_NUM_RADAII,0);
        psi4L2I.resize(BSSN_GW_NUM_RADAII,0);
        psi4L2R_g.resize(BSSN_GW_NUM_RADAII,0);
        psi4L2I_g.resize(BSSN_GW_NUM_RADAII,0);

        if(mesh->isActive())
        {
            const unsigned int rankActive=mesh->getMPIRank();
            const unsigned int npesActive=mesh->getMPICommSize();

            const unsigned int numPts=LEBEDEV_NUM_PTS;

            std::vector<double> domain_coords;
            domain_coords.resize(3*LEBEDEV_NUM_PTS);

            std::vector<double> psi4_real;
            psi4_real.resize(LEBEDEV_NUM_PTS);

            std::vector<double> psi4_imag;
            psi4_imag.resize(LEBEDEV_NUM_PTS);

            std::vector<unsigned int > validIndex;

            Point grid_limits[2];
            Point domain_limits[2];

            grid_limits[0] = Point(bssn::BSSN_OCTREE_MIN[0], bssn::BSSN_OCTREE_MIN[1], bssn::BSSN_OCTREE_MIN[2]);
            grid_limits[1] = Point(bssn::BSSN_OCTREE_MAX[0], bssn::BSSN_OCTREE_MAX[1], bssn::BSSN_OCTREE_MAX[2]);

            domain_limits[0] = Point(bssn::BSSN_COMPD_MIN[0], bssn::BSSN_COMPD_MIN[1], bssn::BSSN_COMPD_MIN[2]);
            domain_limits[1] = Point(bssn::BSSN_COMPD_MAX[0], bssn::BSSN_COMPD_MAX[1], bssn::BSSN_COMPD_MAX[2]);

            for(unsigned int k=0;k<BSSN_GW_NUM_RADAII;k++)
            {

                for(unsigned int pts=0;pts<numPts;pts++)
                {
                    domain_coords[3*pts + 0]=(BSSN_GW_RADAII[k]*sin(LEBEDEV_THETA[pts])*cos(LEBEDEV_PHI[pts]));
                    domain_coords[3*pts + 1]=(BSSN_GW_RADAII[k]*sin(LEBEDEV_THETA[pts])*sin(LEBEDEV_PHI[pts]));
                    domain_coords[3*pts + 2]=(BSSN_GW_RADAII[k]*cos(LEBEDEV_THETA[pts]));
        
                }

                validIndex.clear();
                ot::da::interpolateToCoords(mesh,cVar[bssn::VAR_CONSTRAINT::C_PSI4_REAL],domain_coords.data(),domain_coords.size(), grid_limits, domain_limits, &(*(psi4_real.begin())),validIndex);

                validIndex.clear();
                ot::da::interpolateToCoords(mesh,cVar[bssn::VAR_CONSTRAINT::C_PSI4_IMG],domain_coords.data(),domain_coords.size(), grid_limits, domain_limits, &(*(psi4_imag.begin())),validIndex);

                for(unsigned int index=0;index<validIndex.size();index++)
                {
                    psi4L2R[k] += (psi4_real[validIndex[index]]*psi4_real[validIndex[index]]);
                    psi4L2I[k] += (psi4_imag[validIndex[index]]*psi4_imag[validIndex[index]]);
                }
                    

                for(unsigned int l=0;l<BSSN_GW_NUM_LMODES;l++)
                {
                    for(unsigned int m=0;m<2*BSSN_GW_L_MODES[l]+1;m++)
                    {
                        swsh_coeff[k*TOTAL_MODES+lmOffset[l]+m]=DendroComplex(0.0,0.0);
                        for(unsigned int index=0;index<validIndex.size();index++)
                        {   
                            DendroComplex psi4 = DendroComplex(psi4_real[validIndex[index]],psi4_imag[validIndex[index]]);
                            swsh_coeff[k*TOTAL_MODES+lmOffset[l]+m]+=( psi4 * std::conj(LEBEDEV_SWSH[lmOffset[l]+m][validIndex[index]]) * LEBEDEV_W[validIndex[index]] );
                        }


                        //Eric and I both agreed that r^2 here is not necessary but the 4*M_PI should be there. 
                        // swsh_coeff[k*TOTAL_MODES+lmOffset[l]+m]*=(4*M_PI*BSSN_GW_RADAII[k]*BSSN_GW_RADAII[k]);
                        swsh_coeff[k*TOTAL_MODES+lmOffset[l]+m]*=(4*M_PI);

                    }
                }

            }
        }

        //par::Mpi_Reduce(swsh_coeff,swsh_coeff_g,(BSSN_GW_NUM_RADAII*TOTAL_MODES),MPI_SUM,0,commActive);
        //for(unsigned int k=0;k<BSSN_GW_NUM_RADAII;k++)
        MPI_Reduce(swsh_coeff,swsh_coeff_g,(BSSN_GW_NUM_RADAII*TOTAL_MODES),MPI_DOUBLE_COMPLEX,MPI_SUM,0,commActive);
        MPI_Reduce(&(*(psi4L2R.begin())),&(*(psi4L2R_g.begin())),BSSN_GW_NUM_RADAII,MPI_DOUBLE,MPI_SUM,0,commActive);
        MPI_Reduce(&(*(psi4L2I.begin())),&(*(psi4L2I_g.begin())),BSSN_GW_NUM_RADAII,MPI_DOUBLE,MPI_SUM,0,commActive);

        if(!rankActive)
        {

            int n=BSSN_GW_NUM_RADAII;


            std::ofstream fileGW;
            char fName[256];
            sprintf(fName,"%s_GW_L2.dat",bssn::BSSN_PROFILE_FILE_PREFIX.c_str());
            fileGW.open (fName,std::ofstream::app);

            if(timestep==0)
			{
            	fileGW<<"TimeStep\t"<<" t\t";
                for(unsigned int r=0;r<BSSN_GW_NUM_RADAII;r++)
					fileGW<<"r="<<BSSN_GW_RADAII[r]<<"\t";
                fileGW<<std::endl;

            }

            fileGW<<timestep<<"\t"<<time<<"\t";
            for(unsigned int r=0;r<BSSN_GW_NUM_RADAII;r++)
            {
                psi4L2R_g[r] = sqrt(psi4L2R_g[r]);
                psi4L2I_g[r] = sqrt(psi4L2I_g[r]);
                fileGW<<"("<< psi4L2R_g[r]<<","<<psi4L2I_g[r]<<")\t";
            }
            fileGW<<std::endl;
            fileGW.close();
                


            for(unsigned int l=0;l<BSSN_GW_NUM_LMODES;l++)
            {
                
                for(unsigned int m=0;m<2*BSSN_GW_L_MODES[l]+1;m++)
                {
					// write data to file
                    std::ofstream fileGW;
            		char fName[256];

            		//std::cout<<"l : "<<BSSN_GW_L_MODES[l]<< " m : "<<(int)(m-BSSN_GW_L_MODES[l])<<std::endl;
            		sprintf(fName,"%s_GW_l%d_m%d.dat",bssn::BSSN_PROFILE_FILE_PREFIX.c_str(),BSSN_GW_L_MODES[l],(int)(m-BSSN_GW_L_MODES[l]));
            		fileGW.open (fName,std::ofstream::app);

					if(timestep==0)
					{
						fileGW<<"TimeStep\t"<<" t\t";
						for(unsigned int r=0;r<BSSN_GW_NUM_RADAII;r++)
							fileGW<<"r"<<r<<"\t";

						fileGW<<std::endl;							
					}

                    fileGW.precision(GW::BSSN_GW_OUTPUT_PRECISION);
                    fileGW<<std::scientific;
					fileGW<<timestep<<"\t"<<time<<"\t";
					for(unsigned int r=0;r<BSSN_GW_NUM_RADAII;r++)
						fileGW<<swsh_coeff_g[r*TOTAL_MODES+lmOffset[l]+m]<<"\t";

					fileGW<<std::endl;
					
					fileGW.close();

                }
				
            }

        }

        delete [] swsh_coeff;
        delete [] swsh_coeff_g;

        return;

    }


} // end of namespace GW



#endif

