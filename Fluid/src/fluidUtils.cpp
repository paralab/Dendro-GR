//
// Created by milinda on 7/26/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief Contains utility functions for FLUID simulation.
*/
//

#include "fluidUtils.h"
#include <waveletRefEl.h>

namespace fluid
{

    void readParamFile(const char * fName,MPI_Comm comm)
    {


        json parFile;
        int rank,npes;
        MPI_Comm_rank(comm,&rank);
        MPI_Comm_size(comm,&npes);

        unsigned int vtu_len;
        unsigned int chp_len;
        unsigned int prf_len;

        if(!rank)
        {
            std::ifstream infile(fName);
            if(!infile) {std::cout<<fName<<" parameter file open failed "<<std::endl;}
            infile>>parFile;

            fluid::FLUID_IO_OUTPUT_FREQ=parFile["FLUID_IO_OUTPUT_FREQ"];
            fluid::FLUID_REMESH_TEST_FREQ=parFile["FLUID_REMESH_TEST_FREQ"];
            fluid::FLUID_CHECKPT_FREQ=parFile["FLUID_CHECKPT_FREQ"];
            fluid::FLUID_IO_OUTPUT_GAP=parFile["FLUID_IO_OUTPUT_GAP"];
            fluid::FLUID_VTU_FILE_PREFIX=parFile["FLUID_VTU_FILE_PREFIX"].get<std::string>();
            fluid::FLUID_CHKPT_FILE_PREFIX=parFile["FLUID_CHKPT_FILE_PREFIX"].get<std::string>();
            fluid::FLUID_PROFILE_FILE_PREFIX=parFile["FLUID_PROFILE_FILE_PREFIX"].get<std::string>();
            fluid::FLUID_RESTORE_SOLVER=parFile["FLUID_RESTORE_SOLVER"];
            fluid::FLUID_ENABLE_BLOCK_ADAPTIVITY=parFile["FLUID_ENABLE_BLOCK_ADAPTIVITY"];
            fluid::FLUID_BLK_MIN_X=parFile["FLUID_BLK_MIN_X"];
            fluid::FLUID_BLK_MIN_Y=parFile["FLUID_BLK_MIN_Y"];
            fluid::FLUID_BLK_MIN_Z=parFile["FLUID_BLK_MIN_Z"];
            fluid::FLUID_BLK_MAX_X=parFile["FLUID_BLK_MAX_X"];
            fluid::FLUID_BLK_MAX_Y=parFile["FLUID_BLK_MAX_Y"];
            fluid::FLUID_BLK_MAX_Z=parFile["FLUID_BLK_MAX_Z"];
            fluid::FLUID_DENDRO_GRAIN_SZ=parFile["FLUID_DENDRO_GRAIN_SZ"];
            fluid::FLUID_ASYNC_COMM_K=parFile["FLUID_ASYNC_COMM_K"];
            fluid::FLUID_DENDRO_AMR_FAC=parFile["FLUID_DENDRO_AMR_FAC"];
            fluid::FLUID_LOAD_IMB_TOL=parFile["FLUID_LOAD_IMB_TOL"];
            fluid::FLUID_RK_TIME_BEGIN=parFile["FLUID_RK_TIME_BEGIN"];
            fluid::FLUID_RK_TIME_END=parFile["FLUID_RK_TIME_END"];
            fluid::FLUID_RK45_TIME_STEP_SIZE=parFile["FLUID_RK45_TIME_STEP_SIZE"];
            fluid::FLUID_RK45_DESIRED_TOL=parFile["FLUID_RK45_DESIRED_TOL"];
            fluid::FLUID_DIM=parFile["FLUID_DIM"];
            fluid::FLUID_MAXDEPTH=parFile["FLUID_MAXDEPTH"];
            fluid::FLUID_GRID_MIN_X=parFile["FLUID_GRID_MIN_X"];
            fluid::FLUID_GRID_MAX_X=parFile["FLUID_GRID_MAX_X"];
            fluid::FLUID_GRID_MIN_Y=parFile["FLUID_GRID_MIN_Y"];
            fluid::FLUID_GRID_MAX_Y=parFile["FLUID_GRID_MAX_Y"];
            fluid::FLUID_GRID_MIN_Z=parFile["FLUID_GRID_MIN_Z"];
            fluid::FLUID_GRID_MAX_Z=parFile["FLUID_GRID_MAX_Z"];
            fluid::KO_DISS_SIGMA=parFile["KO_DISS_SIGMA"];

            fluid::FLUID_ID_TYPE=parFile["FLUID_ID_TYPE"];
            fluid::FLUID_RECON_METHOD=parFile["FLUID_RECON_METHOD"];
            if(parFile.find("FLUID_CYL_GAUSSIAN_AXIS")!=parFile.end()){
              fluid::FLUID_CYL_GAUSSIAN_AXIS=parFile["FLUID_CYL_GAUSSIAN_AXIS"];
            }

            fluid::FLUID_WAVELET_TOL=parFile["FLUID_WAVELET_TOL"];

            fluid::FLUID_NUM_REFINE_VARS=parFile["FLUID_NUM_REFINE_VARS"];
            for(unsigned int i=0;i<fluid::FLUID_NUM_REFINE_VARS;i++)
                fluid::FLUID_REFINE_VARIABLE_INDICES[i]=parFile["FLUID_REFINE_VARIABLE_INDICES"][i];

            fluid::FLUID_NUM_EVOL_VARS_VTU_OUTPUT=parFile["FLUID_NUM_EVOL_VARS_VTU_OUTPUT"];

            for(unsigned int i=0;i<fluid::FLUID_NUM_EVOL_VARS_VTU_OUTPUT;i++)
                fluid::FLUID_VTU_OUTPUT_EVOL_INDICES[i]=parFile["FLUID_VTU_OUTPUT_EVOL_INDICES"][i];

            if(parFile.find("FLUID_VTU_Z_SLICE_ONLY")!=parFile.end()){
              fluid::FLUID_VTU_Z_SLICE_ONLY = parFile["FLUID_VTU_Z_SLICE_ONLY"];
            }
            if(parFile.find("FLUID_GAMMA")!=parFile.end()){
              fluid::FLUID_GAMMA = parFile["FLUID_GAMMA"];
            }
            if(parFile.find("FLUID_VACUUM")!=parFile.end()){
              fluid::FLUID_VACUUM = parFile["FLUID_VACUUM"];
              fluid::FLUID_VACUUM_RESET = fluid::FLUID_VACUUM;
            }
            if(parFile.find("FLUID_VACUUM_D")!=parFile.end()){
              fluid::FLUID_VACUUM_D = parFile["FLUID_VACUUM_D"];
              fluid::FLUID_VACUUM_D_RESET = fluid::FLUID_VACUUM_D;
            }
            if(parFile.find("FLUID_VACUUM_TAU")!=parFile.end()){
              fluid::FLUID_VACUUM_TAU = parFile["FLUID_VACUUM_TAU"];
              fluid::FLUID_VACUUM_TAU_RESET = fluid::FLUID_VACUUM_TAU;
            }

            // Duffell-MacFadyen parameters
            if(parFile.find("FLUID_DUFF_MAC_RHO0")!=parFile.end()){
              fluid::FLUID_DUFF_MAC_RHO0 = parFile["FLUID_DUFF_MAC_RHO0"];
            }
            if(parFile.find("FLUID_DUFF_MAC_K0")!=parFile.end()){
              fluid::FLUID_DUFF_MAC_K0 = parFile["FLUID_DUFF_MAC_K0"];
            }
            if(parFile.find("FLUID_DUFF_MAC_K")!=parFile.end()){
              fluid::FLUID_DUFF_MAC_K = parFile["FLUID_DUFF_MAC_K"];
            }
            if(parFile.find("FLUID_DUFF_MAC_RMIN")!=parFile.end()){
              fluid::FLUID_DUFF_MAC_RMIN = parFile["FLUID_DUFF_MAC_RMIN"];
            }
            if(parFile.find("FLUID_DUFF_MAC_R0")!=parFile.end()){
              fluid::FLUID_DUFF_MAC_R0 = parFile["FLUID_DUFF_MAC_R0"];
            }
            if(parFile.find("FLUID_DUFF_MAC_REXP")!=parFile.end()){
              fluid::FLUID_DUFF_MAC_REXP = parFile["FLUID_DUFF_MAC_REXP"];
            }
            if(parFile.find("FLUID_DUFF_MAC_E0")!=parFile.end()){
              fluid::FLUID_DUFF_MAC_E0 = parFile["FLUID_DUFF_MAC_E0"];
            }


            vtu_len=FLUID_VTU_FILE_PREFIX.size();
            chp_len=FLUID_CHKPT_FILE_PREFIX.size();
            prf_len=FLUID_PROFILE_FILE_PREFIX.size();

        }


        par::Mpi_Bcast(&FLUID_IO_OUTPUT_FREQ,1,0,comm);
        par::Mpi_Bcast(&FLUID_REMESH_TEST_FREQ,1,0,comm);
        par::Mpi_Bcast(&FLUID_CHECKPT_FREQ,1,0,comm);
        par::Mpi_Bcast(&FLUID_IO_OUTPUT_GAP,1,0,comm);

        par::Mpi_Bcast(&vtu_len,1,0,comm);
        par::Mpi_Bcast(&chp_len,1,0,comm);
        par::Mpi_Bcast(&prf_len,1,0,comm);

        par::Mpi_Bcast(&FLUID_DENDRO_GRAIN_SZ,1,0,comm);
        par::Mpi_Bcast(&FLUID_DENDRO_AMR_FAC,1,0,comm);
        par::Mpi_Bcast(&FLUID_ASYNC_COMM_K,1,0,comm);

        par::Mpi_Bcast(&FLUID_GAMMA,1,0,comm);
        par::Mpi_Bcast(&FLUID_VTU_Z_SLICE_ONLY,1,0,comm);
        par::Mpi_Bcast(&FLUID_CYL_GAUSSIAN_AXIS,1,0,comm);
        par::Mpi_Bcast(&FLUID_VACUUM,1,0,comm);
        par::Mpi_Bcast(&FLUID_VACUUM_RESET,1,0,comm);
        par::Mpi_Bcast(&FLUID_VACUUM_D,1,0,comm);
        par::Mpi_Bcast(&FLUID_VACUUM_D_RESET,1,0,comm);
        par::Mpi_Bcast(&FLUID_VACUUM_TAU,1,0,comm);
        par::Mpi_Bcast(&FLUID_VACUUM_TAU_RESET,1,0,comm);

        par::Mpi_Bcast(&FLUID_DUFF_MAC_RHO0,1,0,comm);
        par::Mpi_Bcast(&FLUID_DUFF_MAC_K0,1,0,comm);
        par::Mpi_Bcast(&FLUID_DUFF_MAC_K,1,0,comm);
        par::Mpi_Bcast(&FLUID_DUFF_MAC_RMIN,1,0,comm);
        par::Mpi_Bcast(&FLUID_DUFF_MAC_R0,1,0,comm);
        par::Mpi_Bcast(&FLUID_DUFF_MAC_REXP,1,0,comm);
        par::Mpi_Bcast(&FLUID_DUFF_MAC_E0,1,0,comm);

        char vtu_name[vtu_len+1];
        char chp_name[chp_len+1];
        char prf_name[prf_len+1];


        if(!rank)
        {
           for(unsigned int k=0;k<vtu_len;k++)
               vtu_name[k]=FLUID_VTU_FILE_PREFIX[k];

            for(unsigned int k=0;k<chp_len;k++)
                chp_name[k]=FLUID_CHKPT_FILE_PREFIX[k];

            for(unsigned int k=0;k<prf_len;k++)
                prf_name[k]=FLUID_PROFILE_FILE_PREFIX[k];

            vtu_name[vtu_len]='\0';
            chp_name[chp_len]='\0';
            prf_name[prf_len]='\0';

        }


        MPI_Bcast(vtu_name,vtu_len+1,MPI_CHAR,0,comm);
        MPI_Bcast(chp_name,chp_len+1,MPI_CHAR,0,comm);
        MPI_Bcast(prf_name,prf_len+1,MPI_CHAR,0,comm);

        FLUID_VTU_FILE_PREFIX=std::string(vtu_name);
        FLUID_CHKPT_FILE_PREFIX=std::string(chp_name);
        FLUID_PROFILE_FILE_PREFIX=std::string(prf_name);


        par::Mpi_Bcast(&FLUID_RESTORE_SOLVER,1,0,comm);
        par::Mpi_Bcast(&FLUID_ENABLE_BLOCK_ADAPTIVITY,1,0,comm);

        par::Mpi_Bcast(&FLUID_WAVELET_TOL,1,0,comm);

        par::Mpi_Bcast(&FLUID_LOAD_IMB_TOL,1,0,comm);
        par::Mpi_Bcast(&FLUID_RK_TIME_BEGIN,1,0,comm);
        par::Mpi_Bcast(&FLUID_RK_TIME_END,1,0,comm);
        par::Mpi_Bcast(&FLUID_RK45_TIME_STEP_SIZE,1,0,comm);
        par::Mpi_Bcast(&FLUID_RK45_DESIRED_TOL,1,0,comm);
        par::Mpi_Bcast(&FLUID_DIM,1,0,comm);
        par::Mpi_Bcast(&FLUID_MAXDEPTH,1,0,comm);

        par::Mpi_Bcast(&FLUID_ID_TYPE,1,0,comm);
        par::Mpi_Bcast(&FLUID_RECON_METHOD,1,0,comm);

        par::Mpi_Bcast(&FLUID_GRID_MIN_X,1,0,comm);
        par::Mpi_Bcast(&FLUID_GRID_MAX_X,1,0,comm);
        par::Mpi_Bcast(&FLUID_GRID_MIN_Y,1,0,comm);
        par::Mpi_Bcast(&FLUID_GRID_MAX_Y,1,0,comm);
        par::Mpi_Bcast(&FLUID_GRID_MIN_Z,1,0,comm);
        par::Mpi_Bcast(&FLUID_GRID_MAX_Z,1,0,comm);

        par::Mpi_Bcast(&FLUID_BLK_MIN_X,1,0,comm);
        par::Mpi_Bcast(&FLUID_BLK_MIN_Y,1,0,comm);
        par::Mpi_Bcast(&FLUID_BLK_MIN_Z,1,0,comm);

        par::Mpi_Bcast(&FLUID_BLK_MAX_X,1,0,comm);
        par::Mpi_Bcast(&FLUID_BLK_MAX_Y,1,0,comm);
        par::Mpi_Bcast(&FLUID_BLK_MAX_Z,1,0,comm);

        FLUID_OCTREE_MAX[0]=(double )(1u<<fluid::FLUID_MAXDEPTH);
        FLUID_OCTREE_MAX[1]=(double )(1u<<fluid::FLUID_MAXDEPTH);
        FLUID_OCTREE_MAX[2]=(double )(1u<<fluid::FLUID_MAXDEPTH);

        FLUID_COMPD_MIN[0]=FLUID_GRID_MIN_X;
        FLUID_COMPD_MIN[1]=FLUID_GRID_MIN_Y;
        FLUID_COMPD_MIN[2]=FLUID_GRID_MIN_Z;

        FLUID_COMPD_MAX[0]=FLUID_GRID_MAX_X;
        FLUID_COMPD_MAX[1]=FLUID_GRID_MAX_Y;
        FLUID_COMPD_MAX[2]=FLUID_GRID_MAX_Z;


        par::Mpi_Bcast(&KO_DISS_SIGMA, 1, 0, comm);

        par::Mpi_Bcast(&FLUID_NUM_REFINE_VARS,1,0,comm);
        par::Mpi_Bcast(&FLUID_NUM_EVOL_VARS_VTU_OUTPUT,1,0,comm);


        if(FLUID_NUM_REFINE_VARS>FLUID_NUM_PRIM_VARS){std::cout<<"Error[parameter file]: Number of refine variables should be less than number of FLUID_NUM_PRIM_VARS"<<std::endl; exit(0);}
        if(FLUID_NUM_EVOL_VARS_VTU_OUTPUT>FLUID_NUM_PRIM_VARS){std::cout<<"Error[parameter file]: Number of evolution VTU variables should be less than number of FLUID_NUM_PRIM_VARS"<<std::endl; exit(0);}

        par::Mpi_Bcast(FLUID_REFINE_VARIABLE_INDICES,FLUID_NUM_PRIM_VARS,0,comm);
        par::Mpi_Bcast(FLUID_VTU_OUTPUT_EVOL_INDICES,FLUID_NUM_PRIM_VARS,0,comm);

    }



    void initData(const DendroScalar xx1,const DendroScalar yy1,const DendroScalar zz1, DendroScalar *var)
    {

        const DendroScalar x=GRIDX_TO_X(xx1);
        const DendroScalar y=GRIDY_TO_Y(yy1);
        const DendroScalar z=GRIDZ_TO_Z(zz1);
        const DendroScalar r=sqrt(x*x + y*y + z*z);

        //std::cout<<"initData: type="<<fluid::FLUID_ID_TYPE<<", "<<x<<", "<<y<<", "<<z<<std::endl;

        if (fluid::FLUID_ID_TYPE == 0) {
          // For the sake of testing, let's just put a Gaussian here for the time being.
          var[PVAR::V_RHO] = FLUID_VACUUM_RESET + 10.0*exp(-r*r);
          var[PVAR::V_VX]  = 0.0;
          var[PVAR::V_VY]  = 0.0;
          var[PVAR::V_VZ]  = 0.0;
          var[PVAR::V_P]   = FLUID_VACUUM_TAU_RESET + 100.0*pow(var[PVAR::V_RHO],FLUID_GAMMA);

          // 4 - velocity vector
          double v4[3];
          double v3[3]={var[PVAR::V_VX],var[PVAR::V_VY],var[PVAR::V_VZ]};

          double alpha=1.0;
          double Beta[3]={0.0,0.0,0.0};
          double gd[3][3];
          gd[0][0]=1.0;
          gd[0][1]=0.0;
          gd[0][2]=0.0;

          gd[1][0]=0.0;
          gd[1][1]=1.0;
          gd[1][2]=0.0;

          gd[2][0]=0.0;
          gd[2][1]=0.0;
          gd[2][2]=1.0;


          fluid_math::Vvec3_to_Vvec4(v3,v4,gd,alpha,Beta);  

          var[PVAR::V_U1]  = v4[0];  
          var[PVAR::V_U2]  = v4[1];  
          var[PVAR::V_U3]  = v4[2];  
          
        }
       	else if (fluid::FLUID_ID_TYPE == 1) {
          // Let's do a shock tube here
          var[PVAR::V_RHO] = 1.0;
          var[PVAR::V_VX] = 0.0;
          var[PVAR::V_VY] = 0.0;
          var[PVAR::V_VZ] = 0.0;
          if(x > -1.0){
            var[PVAR::V_P] = 0.01;
          }
          else{
            var[PVAR::V_P] = 1000.0;
          }

          double v4[3];
          double v3[3]={var[PVAR::V_VX],var[PVAR::V_VY],var[PVAR::V_VZ]};

          double alpha=1.0;
          double Beta[3]={0.0,0.0,0.0};
          double gd[3][3];
          gd[0][0]=1.0;
          gd[0][1]=0.0;
          gd[0][2]=0.0;

          gd[1][0]=0.0;
          gd[1][1]=1.0;
          gd[1][2]=0.0;

          gd[2][0]=0.0;
          gd[2][1]=0.0;
          gd[2][2]=1.0;


          fluid_math::Vvec3_to_Vvec4(v3,v4,gd,alpha,Beta);  

          var[PVAR::V_U1]  = v4[0];  
          var[PVAR::V_U2]  = v4[1];  
          var[PVAR::V_U3]  = v4[2];  
        }
        else if (fluid::FLUID_ID_TYPE == 2){
          // Let's do a shock tube here
          var[PVAR::V_RHO] = 1.0;
          var[PVAR::V_VX] = 0.0;
          var[PVAR::V_VY] = 0.0;
          var[PVAR::V_VZ] = 0.0;
          if(y > -1.0){
            var[PVAR::V_P] = 0.01;
          }
          else{
            var[PVAR::V_P] = 1000.0;
          }

          double v4[3];
          double v3[3]={var[PVAR::V_VX],var[PVAR::V_VY],var[PVAR::V_VZ]};

          double alpha=1.0;
          double Beta[3]={0.0,0.0,0.0};
          double gd[3][3];
          gd[0][0]=1.0;
          gd[0][1]=0.0;
          gd[0][2]=0.0;

          gd[1][0]=0.0;
          gd[1][1]=1.0;
          gd[1][2]=0.0;

          gd[2][0]=0.0;
          gd[2][1]=0.0;
          gd[2][2]=1.0;


          fluid_math::Vvec3_to_Vvec4(v3,v4,gd,alpha,Beta);  

          var[PVAR::V_U1]  = v4[0];  
          var[PVAR::V_U2]  = v4[1];  
          var[PVAR::V_U3]  = v4[2];  
 
        }
        else if (fluid::FLUID_ID_TYPE == 3){
          // Let's do a shock tube here
          var[PVAR::V_RHO] = 1.0;
          var[PVAR::V_VX] = 0.0;
          var[PVAR::V_VY] = 0.0;
          var[PVAR::V_VZ] = 0.0;
          if(z > -1.0){
            var[PVAR::V_P] = 0.01;
          }
          else{
            var[PVAR::V_P] = 1000.0;
          }

          double v4[3];
          double v3[3]={var[PVAR::V_VX],var[PVAR::V_VY],var[PVAR::V_VZ]};

          double alpha=1.0;
          double Beta[3]={0.0,0.0,0.0};
          double gd[3][3];
          gd[0][0]=1.0;
          gd[0][1]=0.0;
          gd[0][2]=0.0;

          gd[1][0]=0.0;
          gd[1][1]=1.0;
          gd[1][2]=0.0;

          gd[2][0]=0.0;
          gd[2][1]=0.0;
          gd[2][2]=1.0;


          fluid_math::Vvec3_to_Vvec4(v3,v4,gd,alpha,Beta);  

          var[PVAR::V_U1]  = v4[0];  
          var[PVAR::V_U2]  = v4[1];  
          var[PVAR::V_U3]  = v4[2];  
 
        }
        else if (fluid::FLUID_ID_TYPE == 4){
          // We'll set up a Blandford-McKee profile.
          double t = 0.10;
          double rho0 = 7400.66;
          double E = 5.59424;
          double k = 2.0;
          double r0 = 0.01;
          double ls = pow((3.0 - k)*E/(4*M_PI*pow(r0,k)*rho0),1.0/(3.0 - k));
          double fac = 3.0/(3.0 - k);
          double A0 = fac*pow(r0,k)*rho0/pow(ls,k);
          double Ak = A0*pow(ls,k)*(3.0 - k)/3.0;
          double ils = 1.0/ls;
          double eta = 1.0e-10;

          double R = 0.10;
          double Ws = sqrt(pow(R,k)*(17.0 - 4.0*k)*E/(8*M_PI*Ak*t*t*t));
          
          double fac2 = (4.0 + k)/(6.0*(4.0 - k));
          double fac3 = 2.0*(4.0 - k)*Ws*Ws;
          double ifac = 1.0/(3.0*sqrt(2.0));
          double isqrt2 = 1.0/sqrt(2.0);

          double s = sqrt(x*x + y*y);
          double sthcph = x/(r+1.0e-15);
          double sthsph = y/(r+1.0e-15);
          double cth = z/(r+1.0e-15);

          // If we are effectively at the shock, then it needs to be represented in our data.
          // FIXME: Adjust to work with grid spacing rather than arbitrary data.
          if(fabs(r-R) < 1e-15){
            double W = Ws;
            var[PVAR::V_RHO] = fac*Ak*pow(R + 1.0e-15,-k)*Ws;
            var[PVAR::V_P] = ifac*var[PVAR::V_RHO]*Ws;
            double vt = sqrt(W*W - 1.0)/W;

            // Set the velocity
            var[PVAR::V_VX] = vt*sthcph;
            var[PVAR::V_VY] = vt*sthsph;
            var[PVAR::V_VZ] = vt*cth;
          }
          else if(r < R){
            double chi = 1.0 + fac3*(1.0 - r/R);
            double W = fmax(1.0,isqrt2*Ws/sqrt(chi));
            double vt = sqrt(W*W - 1.0)/W;
            // Set the density
            var[PVAR::V_RHO] = pow(2.0,1.5)*Ak*pow(R,-k)*Ws*pow(chi,-(10.0 - 3.0*k)/(2.0*(4.0 - k)));
            // Set the pressure
            var[PVAR::V_P] = ifac*var[PVAR::V_RHO]*Ws*pow(chi,-fac2);
            // Set the velocity
            var[PVAR::V_VX] = vt*sthcph;
            var[PVAR::V_VY] = vt*sthsph;
            var[PVAR::V_VZ] = vt*cth;
          }
          else{
            var[PVAR::V_RHO] = Ak*pow(r,-k);
            var[PVAR::V_P] = var[PVAR::V_RHO]*eta;
            var[PVAR::V_VX] = 0.0;
            var[PVAR::V_VY] = 0.0;
            var[PVAR::V_VZ] = 0.0;
          }
          double v4[3];
          double v3[3]={var[PVAR::V_VX],var[PVAR::V_VY],var[PVAR::V_VZ]};

          double alpha=1.0;
          double Beta[3]={0.0,0.0,0.0};
          double gd[3][3];
          gd[0][0]=1.0;
          gd[0][1]=0.0;
          gd[0][2]=0.0;

          gd[1][0]=0.0;
          gd[1][1]=1.0;
          gd[1][2]=0.0;

          gd[2][0]=0.0;
          gd[2][1]=0.0;
          gd[2][2]=1.0;


          fluid_math::Vvec3_to_Vvec4(v3,v4,gd,alpha,Beta);  

          var[PVAR::V_U1]  = v4[0];  
          var[PVAR::V_U2]  = v4[1];  
          var[PVAR::V_U3]  = v4[2];  
        }
        else if (fluid::FLUID_ID_TYPE == 5){
          // Constant initial data.
          var[PVAR::V_RHO] = 1;
          var[PVAR::V_VX ] = 0.0;
          var[PVAR::V_VY ] = 0.0;
          var[PVAR::V_VZ ] = 0.0;
          var[PVAR::V_P  ] = 0.01;
          var[PVAR::V_U1 ] = 0.0;
          var[PVAR::V_U2 ] = 0.0;
          var[PVAR::V_U3 ] = 0.0;
        }
        else if (fluid::FLUID_ID_TYPE == 6){
          // A cylindrical Gaussian.
          DendroScalar s;
          switch(FLUID_CYL_GAUSSIAN_AXIS){
            case 0:
              s = std::sqrt(z*z + y*y);
              break;
            case 1:
              s = std::sqrt(x*x + z*z);
              break;
            case 2:
              s = std::sqrt(x*x + y*y);
              break;
            default:
              s = std::sqrt(x*x + y*y);
          }
          var[PVAR::V_RHO] = FLUID_VACUUM_RESET + 10.0*exp(-s*s);
          var[PVAR::V_VX]  = 0.0;
          var[PVAR::V_VY]  = 0.0;
          var[PVAR::V_VZ]  = 0.0;
          var[PVAR::V_P]   = FLUID_VACUUM_TAU_RESET + 100.0*pow(var[PVAR::V_RHO],FLUID_GAMMA);

          // 4 - velocity vector
          double v4[3];
          double v3[3]={var[PVAR::V_VX],var[PVAR::V_VY],var[PVAR::V_VZ]};

          double alpha=1.0;
          double Beta[3]={0.0,0.0,0.0};
          double gd[3][3];
          gd[0][0]=1.0;
          gd[0][1]=0.0;
          gd[0][2]=0.0;

          gd[1][0]=0.0;
          gd[1][1]=1.0;
          gd[1][2]=0.0;

          gd[2][0]=0.0;
          gd[2][1]=0.0;
          gd[2][2]=1.0;


          fluid_math::Vvec3_to_Vvec4(v3,v4,gd,alpha,Beta);  

          var[PVAR::V_U1]  = v4[0];  
          var[PVAR::V_U2]  = v4[1];  
          var[PVAR::V_U3]  = v4[2];  

        }
        else if (fluid::FLUID_ID_TYPE == 7){
          // Parameters, generalize later.
          DendroScalar rho0 = FLUID_DUFF_MAC_RHO0;
          DendroScalar k0 = FLUID_DUFF_MAC_K0;
          DendroScalar k = FLUID_DUFF_MAC_K;
          DendroScalar rmin = FLUID_DUFF_MAC_RMIN;
          DendroScalar r0 = FLUID_DUFF_MAC_R0;
          DendroScalar rexp = FLUID_DUFF_MAC_REXP;
          DendroScalar e0 = FLUID_DUFF_MAC_E0;
          
          DendroScalar r = std::sqrt(x*x + y*y + z*z);

          // Set the density profile.
          if(r < rmin){
            var[PVAR::V_RHO] = FLUID_VACUUM_RESET + rho0*std::pow(r0/rmin,k0);
          }
          else if(r < r0){
            var[PVAR::V_RHO] = FLUID_VACUUM_RESET + rho0*std::pow(r0/r,k0);
          }
          else{
            var[PVAR::V_RHO] = FLUID_VACUUM_RESET + rho0*std::pow(r0/r,k);
          }

          // Set the pressure profile.
          if(r < rexp){
            var[PVAR::V_P] = FLUID_VACUUM_TAU_RESET + var[PVAR::V_RHO]*e0/3.0;
          }
          else{
            var[PVAR::V_P] = FLUID_VACUUM_TAU_RESET + var[PVAR::V_RHO]*1e-6;
          }

          var[PVAR::V_VX] = 0.0;
          var[PVAR::V_VY] = 0.0;
          var[PVAR::V_VZ] = 0.0;

          var[PVAR::V_U1] = 0.0;
          var[PVAR::V_U2] = 0.0;
          var[PVAR::V_U3] = 0.0;
        }
        else if(fluid::FLUID_ID_TYPE == 8){
          // Parameters, generalize later.
          DendroScalar rho0 = FLUID_DUFF_MAC_RHO0;
          DendroScalar k0 = FLUID_DUFF_MAC_K0;
          DendroScalar k = FLUID_DUFF_MAC_K;
          DendroScalar rmin = FLUID_DUFF_MAC_RMIN;
          DendroScalar r0 = FLUID_DUFF_MAC_R0;
          DendroScalar rexp = FLUID_DUFF_MAC_REXP;
          DendroScalar e0 = FLUID_DUFF_MAC_E0;

          // Shell parameters, generalize later.
          DendroScalar A = 1.0e3;
          DendroScalar rshell_min = 0.2;
          DendroScalar rshell_max = 0.2015;
          unsigned int angle_dependence = 0;

          DendroScalar r = std::sqrt(x*x + y*y + z*z);
          DendroScalar s = std::sqrt(x*x + y*y);

          DendroScalar rho;

          // Set the density profile.
          if(r < rmin){
            var[PVAR::V_RHO] = FLUID_VACUUM_RESET + rho0*std::pow(r0/rmin,k0);
            rho = var[PVAR::V_RHO];
          }
          else if (r < r0){
            var[PVAR::V_RHO] = FLUID_VACUUM_RESET + rho0*std::pow(r0/r,k0);
            rho = var[PVAR::V_RHO];
          }
          else if (rshell_min <= r && r <= rshell_max){
            double Ashell;
            
            // Angular dependence of the shell.
            switch(angle_dependence){
              case 0:
                // Shell has no angular dependence; spherically symmetric.
                Ashell = A;
                break;
              case 1:
                // Shell is dipolar and open along the poles.
                Ashell = 1.0 + (A - 1.0)*(s/r)*(s/r);
                break;
              case 2:
                // Shell is dipolar and open along the equator.
                Ashell = 1.0 + (A - 1.0)*(z/r)*(z/r);
                break;
              case 3:
                // Shell is dipolar and open along the poles with an additional
                // octopole angle dependence around the equator.
                double cosphi = x/s;
                Ashell = 1.0 + (A - 1.0)*(s/r)*(s/r)*pow(2.5*cosphi*cosphi*cosphi - 1.5*cosphi, 2);
                break;
            }

            var[PVAR::V_RHO] = FLUID_VACUUM_RESET + Ashell*rho0*std::pow(r0/r, k);
            rho = FLUID_VACUUM_RESET + rho0*std::pow(r0/r, k);
          }
          else{
            var[PVAR::V_RHO] = FLUID_VACUUM_RESET + rho0*std::pow(r0/r,k);
            rho = var[PVAR::V_RHO];
          }

          // Set the pressure profile.
          if(r < rexp){
            var[PVAR::V_P] = FLUID_VACUUM_TAU_RESET + rho*e0/3.0;
          }
          else{
            var[PVAR::V_P] = FLUID_VACUUM_TAU_RESET + rho*1e-6;
          }

          var[PVAR::V_VX] = 0.0;
          var[PVAR::V_VY] = 0.0;
          var[PVAR::V_VZ] = 0.0;

          var[PVAR::V_U1] = 0.0;
          var[PVAR::V_U2] = 0.0;
          var[PVAR::V_U3] = 0.0;
        }


    }

    void solToroidalDipole(const double t, const double xx1, const double yy1, const double zz1, double *var)
    {
      // This is a relic from copying the Maxwell files. It doesn't really
      // have a place here, but I don't want to remove it until I have a
      // chance to look at what the code is really doing that calls it.

    }

    void blockAdaptiveOctree(std::vector<ot::TreeNode>& tmpNodes,const Point& pt_min,const Point & pt_max,const unsigned int regLev,const unsigned int maxDepth,MPI_Comm comm)
    {
        int rank,npes;
        MPI_Comm_size(comm,&npes);
        MPI_Comm_rank(comm,&rank);

        double pt_g_min[3];
        double pt_g_max[3];

        pt_g_min[0]=X_TO_GRIDX(pt_min.x());
        pt_g_min[1]=Y_TO_GRIDY(pt_min.y());
        pt_g_min[2]=Z_TO_GRIDZ(pt_min.z());

        pt_g_max[0]=X_TO_GRIDX(pt_max.x());
        pt_g_max[1]=Y_TO_GRIDY(pt_max.y());
        pt_g_max[2]=Z_TO_GRIDZ(pt_max.z());

        assert(pt_g_min[0]>=0 && pt_g_min[0]<(1u<<maxDepth));
        assert(pt_g_min[1]>=0 && pt_g_min[1]<(1u<<maxDepth));
        assert(pt_g_min[2]>=0 && pt_g_min[2]<(1u<<maxDepth));

        assert(pt_g_max[0]>=0 && pt_g_max[0]<(1u<<maxDepth));
        assert(pt_g_max[1]>=0 && pt_g_max[1]<(1u<<maxDepth));
        assert(pt_g_max[2]>=0 && pt_g_max[2]<(1u<<maxDepth));


        unsigned int xRange_b,xRange_e;
        unsigned int yRange_b=pt_g_min[1],yRange_e=pt_g_max[1];
        unsigned int zRange_b=pt_g_min[2],zRange_e=pt_g_max[2];

        xRange_b=pt_g_min[0];//(rank*(pt_g_max[0]-pt_g_min[0]))/npes + pt_g_min[0];
        xRange_e=pt_g_max[1];//((rank+1)*(pt_g_max[0]-pt_g_min[0]))/npes + pt_g_min[0];

        unsigned int stepSz=1u<<(maxDepth-regLev);

       /* std::cout<<" x min: "<<xRange_b<<" x_max: "<<xRange_e<<std::endl;
        std::cout<<" y min: "<<yRange_b<<" y_max: "<<yRange_e<<std::endl;
        std::cout<<" z min: "<<zRange_b<<" z_max: "<<zRange_e<<std::endl;*/


        for(unsigned int x=xRange_b;x<xRange_e;x+=stepSz)
            for(unsigned int y=yRange_b;y<yRange_e;y+=stepSz)
              for(unsigned int z=zRange_b;z<zRange_e;z+=stepSz)
                  tmpNodes.push_back(ot::TreeNode(x,y,z,regLev,m_uiDim,maxDepth));


        return ;


    }

    double computeWTol(double x,double y,double z,double tolMin)
    {
       return fluid::FLUID_WAVELET_TOL;
    }

}// end of namespace fluid





namespace fluid
{

    namespace timer
    {
        void initFlops()
        {
            total_runtime.start();
            t_f2o.start();
            t_cons.start();
            t_bal.start();
            t_mesh.start();
            t_rkSolve.start();
            t_ghostEx_sync.start();
            t_unzip_sync.start();

            for(unsigned int i=0;i<NUM_FACES;i++)
                dendro::timer::t_unzip_sync_face[i].start();

            dendro::timer::t_unzip_async_internal.start();
            dendro::timer::t_unzip_sync_edge.start();
            dendro::timer::t_unzip_sync_vtex.start();
            dendro::timer::t_unzip_p2c.start();
            dendro::timer::t_unzip_sync_nodalval.start();
            dendro::timer::t_unzip_sync_cpy.start();
            dendro::timer::t_unzip_sync_f_c1.start();
            dendro::timer::t_unzip_sync_f_c2.start();
            dendro::timer::t_unzip_sync_f_c3.start();

            t_unzip_async.start();
            dendro::timer::t_unzip_async_comm.start();

            dendro::timer::t_unzip_async_internal.start();
            dendro::timer::t_unzip_async_external.start();
            dendro::timer::t_unzip_async_comm.start();
            t_deriv.start();
            t_rhs.start();

            t_bdyc.start();

            t_zip.start();
            t_rkStep.start();
            t_isReMesh.start();
            t_gridTransfer.start();
            t_ioVtu.start();
            t_ioCheckPoint.start();
        }

        void resetSnapshot()
        {

            total_runtime.snapreset();
            t_f2o.snapreset();
            t_cons.snapreset();
            t_bal.snapreset();
            t_mesh.snapreset();
            t_rkSolve.snapreset();
            t_ghostEx_sync.snapreset();
            t_unzip_sync.snapreset();

            for(unsigned int i=0;i<NUM_FACES;i++)
                dendro::timer::t_unzip_sync_face[i].snapreset();

            dendro::timer::t_unzip_sync_internal.snapreset();
            dendro::timer::t_unzip_sync_edge.snapreset();
            dendro::timer::t_unzip_sync_vtex.snapreset();
            dendro::timer::t_unzip_p2c.snapreset();
            dendro::timer::t_unzip_sync_nodalval.snapreset();
            dendro::timer::t_unzip_sync_cpy.snapreset();

            dendro::timer::t_unzip_sync_f_c1.snapreset();
            dendro::timer::t_unzip_sync_f_c2.snapreset();
            dendro::timer::t_unzip_sync_f_c3.snapreset();

            t_unzip_async.snapreset();
            dendro::timer::t_unzip_async_internal.snapreset();
            dendro::timer::t_unzip_async_external.snapreset();
            dendro::timer::t_unzip_async_comm.snapreset();

            t_deriv.snapreset();
            t_rhs.snapreset();

            t_bdyc.snapreset();

            t_zip.snapreset();
            t_rkStep.snapreset();
            t_isReMesh.snapreset();
            t_gridTransfer.snapreset();
            t_ioVtu.snapreset();
            t_ioCheckPoint.snapreset();

        }

        void profileInfo(const char* filePrefix,const ot::Mesh* pMesh)
        {


            int activeRank,activeNpes,globalRank,globalNpes;

            MPI_Comm commActive;
            MPI_Comm commGlobal;

            if(pMesh->isActive())
            {
                commActive=pMesh->getMPICommunicator();
                activeRank=pMesh->getMPIRank();
                activeNpes=pMesh->getMPICommSize();
            }

            globalRank=pMesh->getMPIRankGlobal();
            globalNpes=pMesh->getMPICommSizeGlobal();
            commGlobal=pMesh->getMPIGlobalCommunicator();

            double t_stat;
            double t_stat_g[3];

            const char separator    = ' ';
            const int nameWidth     = 30;
            const int numWidth      = 10;

            char fName[256];
            std::ofstream outfile;

            DendroIntL localSz,globalSz;

            if(!activeRank)
            {
                sprintf(fName,"%s_final.prof",filePrefix);
                outfile.open (fName);
                if(outfile.fail()) {std::cout<<fName<<" file open failed "<<std::endl; return ;}

                outfile<<"active npes : "<<activeNpes<<std::endl;
                outfile<<"global npes : "<<globalNpes<<std::endl;
                outfile<<"partition tol : "<<fluid::FLUID_LOAD_IMB_TOL<<std::endl;
                outfile<<"wavelet tol : "<<fluid::FLUID_WAVELET_TOL<<std::endl;
                outfile<<"maxdepth : "<<fluid::FLUID_MAXDEPTH<<std::endl;

            }

            MPI_Comm comm=commActive;
            unsigned int rank =activeRank;

            localSz=pMesh->getNumLocalMeshElements();
            par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);
            if(!rank)outfile<<"Elements : "<<globalSz<<std::endl;

            localSz=pMesh->getNumLocalMeshNodes();
            par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);
            if(!rank)outfile<<"DOG(zip) : "<<globalSz<<std::endl;

            localSz=pMesh->getDegOfFreedomUnZip();
            par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);
            if(!rank)outfile<<"DOG(unzip) : "<<globalSz<<std::endl;


            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) << "step";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) << "min(s)";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) << "mean(s)";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) << "max(s)"<<std::endl;


            t_stat=total_runtime.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"+runtime(s)";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_f2o.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<" ++f2o";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_cons.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<" ++construction";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_rkSolve.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<" ++rkSolve";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_bal.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --2:1 balance";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_mesh.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --mesh";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_rkStep.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --rkstep";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_ghostEx_sync.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --ghostExchge.";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_unzip_sync.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_sync";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_unzip_async.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  ++unzip_async";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

#ifdef ENABLE_DENDRO_PROFILE_COUNTERS
            t_stat=dendro::timer::t_unzip_async_internal.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_internal";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=dendro::timer::t_unzip_async_external.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_external";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=dendro::timer::t_unzip_async_comm.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_comm (comm) ";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;
#endif

            t_stat=t_deriv.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --deriv ";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=t_rhs.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --compute_rhs ";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=t_bdyc.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --boundary con ";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;



            t_stat=t_zip.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --zip";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_ioVtu.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --vtu";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=t_ioCheckPoint.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --checkpoint";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            if(!rank) outfile.close();




        }


        void profileInfoIntermediate(const char* filePrefix,const ot::Mesh* pMesh,const unsigned int currentStep)
        {

            int activeRank,activeNpes,globalRank,globalNpes;

            MPI_Comm commActive;
            MPI_Comm commGlobal;

            if(pMesh->isActive())
            {
                commActive=pMesh->getMPICommunicator();
                activeRank=pMesh->getMPIRank();
                activeNpes=pMesh->getMPICommSize();
            }

            globalRank=pMesh->getMPIRankGlobal();
            globalNpes=pMesh->getMPICommSizeGlobal();
            commGlobal=pMesh->getMPIGlobalCommunicator();

            double t_stat;
            double t_stat_g[3];

            const char separator    = ' ';
            const int nameWidth     = 30;
            const int numWidth      = 10;

            char fName[256];
            std::ofstream outfile;

            DendroIntL localSz,globalSz;

            DendroIntL ghostElements;
            DendroIntL localElements;

            DendroIntL ghostNodes;
            DendroIntL localNodes;

            DendroIntL totalSendNode;
            DendroIntL totalRecvNode;

            DendroIntL numCalls;


#ifdef FLUID_PROFILE_HUMAN_READABLE
            if(!activeRank)
            {
                sprintf(fName,"%s_im.prof",filePrefix);
                outfile.open (fName,std::fstream::app);
                if(outfile.fail()) {std::cout<<fName<<" file open failed "<<std::endl; return ;}

                outfile<<"active npes : "<<activeNpes<<std::endl;
                outfile<<"global npes : "<<globalNpes<<std::endl;
                outfile<<"current step : "<<currentStep<<std::endl;
                outfile<<"partition tol : "<<fluid::FLUID_LOAD_IMB_TOL<<std::endl;
                outfile<<"wavelet tol : "<<fluid::FLUID_WAVELET_TOL<<std::endl;
                outfile<<"maxdepth : "<<fluid::FLUID_MAXDEPTH<<std::endl;

            }

            MPI_Comm comm=commActive;
            unsigned int rank =activeRank;

            localSz=pMesh->getNumLocalMeshElements();
            par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);
            if(!rank)outfile<<"Elements : "<<globalSz<<std::endl;

            localSz=pMesh->getNumLocalMeshNodes();
            par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);
            if(!rank)outfile<<"DOG(zip) : "<<globalSz<<std::endl;

            localSz=pMesh->getDegOfFreedomUnZip();
            par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);
            if(!rank)outfile<<"DOG(unzip) : "<<globalSz<<std::endl;


            ghostElements=pMesh->getNumPreGhostElements()+pMesh->getNumPostGhostElements();
            localElements=pMesh->getNumLocalMeshElements();

            ghostNodes=pMesh->getNumPreMeshNodes()+pMesh->getNumPostMeshNodes();
            localNodes=pMesh->getNumLocalMeshNodes();

            if(!rank)outfile<<"========================= MESH ======================================================================= "<<std::endl;

            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) << "step";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) << "min(#)";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) << "mean(#)";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) << "max(#)"<<std::endl;

            t_stat=ghostElements;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"ghost Elements";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=localElements;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"local Elements";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=ghostNodes;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"ghost Nodes";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=localNodes;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"local Nodes";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=pMesh->getGhostExcgTotalSendNodeCount();
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"send Nodes";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=pMesh->getGhostExcgTotalRecvNodeCount();
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"recv Nodes";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            if(!rank)outfile<<"========================= RUNTIME =================================================================== "<<std::endl;
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) << "step";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) << "min(s)";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) << "mean(s)";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) << "max(s)"<<std::endl;




            /* t_stat=total_runtime.seconds;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"+runtime(s)";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_f2o.seconds;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<" ++f2o";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_cons.seconds;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<" ++construction";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_rkSolve.seconds;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<" ++rkSolve";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;*/


            t_stat=t_bal.snap;
            //numCalls=t_bal.num_calls;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  ++2:1 balance";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_mesh.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  ++mesh";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_rkStep.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  ++rkstep";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_ghostEx_sync.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  ++ghostExchge.";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_unzip_sync.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  ++unzip_sync";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=t_unzip_async.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  ++unzip_async";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;



            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS

            t_stat=dendro::timer::t_unzip_async_comm.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_comm_wait (comm) ";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=dendro::timer::t_unzip_sync_nodalval.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_nodalVal";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=dendro::timer::t_unzip_sync_f_c1.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --t_unzip_sync_f_c1";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=dendro::timer::t_unzip_sync_f_c2.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --t_unzip_sync_f_c2";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=dendro::timer::t_unzip_sync_f_c3.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --t_unzip_sync_f_c3";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=dendro::timer::t_unzip_sync_cpy.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --t_unzip_sync_cpy";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=dendro::timer::t_unzip_sync_internal.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_sync_internal";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=dendro::timer::t_unzip_sync_face[0].snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_sync_face_left";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=dendro::timer::t_unzip_sync_face[1].snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_sync_face_right";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=dendro::timer::t_unzip_sync_face[2].snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_sync_face_down";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=dendro::timer::t_unzip_sync_face[3].snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_sync_face_up";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=dendro::timer::t_unzip_sync_face[4].snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_sync_face_back";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=dendro::timer::t_unzip_sync_face[5].snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_sync_face_front";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=dendro::timer::t_unzip_sync_edge.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_sync_edge";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=dendro::timer::t_unzip_sync_vtex.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_sync_vtex";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=dendro::timer::t_unzip_p2c.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_p2c";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;
            #endif

            /*
            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
            t_stat=t_unzip_async_internal.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_internal";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=t_unzip_async_external.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_external";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_unzip_async_comm.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_comm (comm) ";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;
            #endif
            */
            t_stat=t_isReMesh.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  ++isReMesh";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=t_gridTransfer.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  ++gridTransfer";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;



            t_stat=t_deriv.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  ++deriv ";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=t_rhs.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  ++compute_rhs ";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=t_bdyc.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  ++boundary con ";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;



            t_stat=t_zip.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  ++zip";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_ioVtu.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  ++vtu";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=t_ioCheckPoint.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  ++checkpoint";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            if(!rank) outfile.close();
#else

            if(!activeRank)
            {
                sprintf(fName,"%s_im.prof",filePrefix);
                outfile.open (fName,std::fstream::app);
                if(outfile.fail()) {std::cout<<fName<<" file open failed "<<std::endl; return ;}

                //writes the header
                if(currentStep==0)
                 outfile<<"step\t act_npes\t glb_npes\t part_tol\t wave_tol\t maxdepth\t numOcts\t dof_zip\t dof_unzip\t"<<\
                 "element_ghost_min\t element_ghost_mean\t element_ghost_max\t"<<\
                 "element_local_min\t element_local_mean\t element_local_max\t"<<\
                 "nodes_local_min\t nodes_local_mean\t nodes_local|max\t"<<\
                 "send_nodes_min\t send_nodes_mean\t send_nodes_max\t"<<\
                 "recv_nodes_min\t recv_nodes_mean\t recv_nodes_max\t"<<\
                 "bal_min\t bal_mean\t bal_max\t"<<\
                 "mesh_min\t mesh_mean\t mesh_max\t"<<\
                 "rkstep_min\t rkstep_mean\t rkstep_max\t"<<\
                 "ghostEx_min\t ghostEx_mean\t ghostEx_max\t"<<\
                 "unzip_sync_min\t unzip_sync_mean\t unzip_sync_max\t"<<\
                 "unzip_async_min\t unzip_async_mean\t unzip_async_max\t"<<\
                 "unzip_async_wait_min\t unzip_async_wait_mean\t unzip_async_wait_max\t"<<\
                 "isRemesh_min\t isRemesh_mean\t isRemesh_max\t"<<\
                 "GT_min\t GT_mean\t GT_max\t"<<\
                 "deriv_min\t deriv_mean\t deriv_max\t"<<\
                 "rhs_min\t rhs_mean\t rhs_max\t"<<std::endl;

            }

            MPI_Comm comm=commActive;
            unsigned int rank =activeRank;

            if(!rank) outfile<<currentStep<<"\t ";
            if(!rank) outfile<<activeNpes<<"\t ";
            if(!rank) outfile<<globalNpes<<"\t ";
            if(!rank) outfile<<fluid::FLUID_LOAD_IMB_TOL<<"\t ";
            if(!rank) outfile<<fluid::FLUID_WAVELET_TOL<<"\t ";
            if(!rank) outfile<<fluid::FLUID_MAXDEPTH<<"\t ";

            localSz=pMesh->getNumLocalMeshElements();
            par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);
            if(!rank)outfile<<globalSz<<"\t ";

            localSz=pMesh->getNumLocalMeshNodes();
            par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);
            if(!rank)outfile<<globalSz<<"\t ";

            localSz=pMesh->getDegOfFreedomUnZip();
            par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);
            if(!rank)outfile<<globalSz<<"\t ";

            ghostElements=pMesh->getNumPreGhostElements()+pMesh->getNumPostGhostElements();
            localElements=pMesh->getNumLocalMeshElements();

            t_stat=ghostElements;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

            t_stat=localElements;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

            ghostNodes=pMesh->getNumPreMeshNodes()+pMesh->getNumPostMeshNodes();
            localNodes=pMesh->getNumLocalMeshNodes();

            /*t_stat=ghostNodes;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";*/

            t_stat=localNodes;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

            t_stat=pMesh->getGhostExcgTotalSendNodeCount();
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

            t_stat=pMesh->getGhostExcgTotalRecvNodeCount();
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

            t_stat=t_bal.snap;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

            t_stat=t_mesh.snap;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

            t_stat=t_rkStep.snap;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

            t_stat=t_ghostEx_sync.snap;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

            t_stat=t_unzip_sync.snap;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

            t_stat=t_unzip_async.snap;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

            t_stat=dendro::timer::t_unzip_async_comm.snap;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

            t_stat=t_isReMesh.snap;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

            t_stat=t_gridTransfer.snap;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

            t_stat=t_deriv.snap;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

            t_stat=t_rhs.snap;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

            if(!rank) outfile<<std::endl;
            if(!rank) outfile.close();
#endif




        }


    }

}
