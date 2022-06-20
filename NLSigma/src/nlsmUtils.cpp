//
// Created by milinda on 7/26/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief Contains utility functions for NLSM simulation.
*/
//

#include "nlsmUtils.h"

namespace nlsm
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

            if(parFile.find("NLSM_ELE_ORDER")!=parFile.end())
                nlsm::NLSM_ELE_ORDER = parFile["NLSM_ELE_ORDER"];

            nlsm::NLSM_IO_OUTPUT_FREQ=parFile["NLSM_IO_OUTPUT_FREQ"];
            nlsm::NLSM_REMESH_TEST_FREQ=parFile["NLSM_REMESH_TEST_FREQ"];
            nlsm::NLSM_CHECKPT_FREQ=parFile["NLSM_CHECKPT_FREQ"];
            nlsm::NLSM_IO_OUTPUT_GAP=parFile["NLSM_IO_OUTPUT_GAP"];
            nlsm::NLSM_VTU_FILE_PREFIX=parFile["NLSM_VTU_FILE_PREFIX"].get<std::string>();
            nlsm::NLSM_CHKPT_FILE_PREFIX=parFile["NLSM_CHKPT_FILE_PREFIX"].get<std::string>();
            nlsm::NLSM_PROFILE_FILE_PREFIX=parFile["NLSM_PROFILE_FILE_PREFIX"].get<std::string>();
            nlsm::NLSM_RESTORE_SOLVER=parFile["NLSM_RESTORE_SOLVER"];
            nlsm::NLSM_ENABLE_BLOCK_ADAPTIVITY=parFile["NLSM_ENABLE_BLOCK_ADAPTIVITY"];
            nlsm::NLSM_BLK_MIN_X=parFile["NLSM_BLK_MIN_X"];
            nlsm::NLSM_BLK_MIN_Y=parFile["NLSM_BLK_MIN_Y"];
            nlsm::NLSM_BLK_MIN_Z=parFile["NLSM_BLK_MIN_Z"];
            nlsm::NLSM_BLK_MAX_X=parFile["NLSM_BLK_MAX_X"];
            nlsm::NLSM_BLK_MAX_Y=parFile["NLSM_BLK_MAX_Y"];
            nlsm::NLSM_BLK_MAX_Z=parFile["NLSM_BLK_MAX_Z"];
            nlsm::NLSM_DENDRO_GRAIN_SZ=parFile["NLSM_DENDRO_GRAIN_SZ"];
            nlsm::NLSM_ASYNC_COMM_K=parFile["NLSM_ASYNC_COMM_K"];
            nlsm::NLSM_DENDRO_AMR_FAC=parFile["NLSM_DENDRO_AMR_FAC"];
            nlsm::NLSM_LOAD_IMB_TOL=parFile["NLSM_LOAD_IMB_TOL"];
            nlsm::NLSM_RK45_TIME_BEGIN=parFile["NLSM_RK45_TIME_BEGIN"];
            nlsm::NLSM_RK45_TIME_END=parFile["NLSM_RK45_TIME_END"];
            nlsm::NLSM_RK45_TIME_STEP_SIZE=parFile["NLSM_RK45_TIME_STEP_SIZE"];
            nlsm::NLSM_RK45_DESIRED_TOL=parFile["NLSM_RK45_DESIRED_TOL"];
            nlsm::NLSM_DIM=parFile["NLSM_DIM"];
            nlsm::NLSM_MAXDEPTH=parFile["NLSM_MAXDEPTH"];
            nlsm::NLSM_GRID_MIN_X=parFile["NLSM_GRID_MIN_X"];
            nlsm::NLSM_GRID_MAX_X=parFile["NLSM_GRID_MAX_X"];
            nlsm::NLSM_GRID_MIN_Y=parFile["NLSM_GRID_MIN_Y"];
            nlsm::NLSM_GRID_MAX_Y=parFile["NLSM_GRID_MAX_Y"];
            nlsm::NLSM_GRID_MIN_Z=parFile["NLSM_GRID_MIN_Z"];
            nlsm::NLSM_GRID_MAX_Z=parFile["NLSM_GRID_MAX_Z"];
            nlsm::KO_DISS_SIGMA=parFile["KO_DISS_SIGMA"];

            nlsm::NLSM_ID_TYPE=parFile["NLSM_ID_TYPE"];
            nlsm::NLSM_ID_AMP1=parFile["NLSM_ID_AMP1"];
            nlsm::NLSM_ID_AMP2=parFile["NLSM_ID_AMP2"];
            nlsm::NLSM_ID_DELTA1=parFile["NLSM_ID_DELTA1"];
            nlsm::NLSM_ID_DELTA2=parFile["NLSM_ID_DELTA2"];
            nlsm::NLSM_ID_XC1=parFile["NLSM_ID_XC1"];
            nlsm::NLSM_ID_YC1=parFile["NLSM_ID_YC1"];
            nlsm::NLSM_ID_ZC1=parFile["NLSM_ID_ZC1"];
            nlsm::NLSM_ID_XC2=parFile["NLSM_ID_XC2"];
            nlsm::NLSM_ID_YC2=parFile["NLSM_ID_YC2"];
            nlsm::NLSM_ID_ZC2=parFile["NLSM_ID_ZC2"];
            
            nlsm::NLSM_ID_EPSX1=parFile["NLSM_ID_EPSX1"];
            nlsm::NLSM_ID_EPSY1=parFile["NLSM_ID_EPSY1"];
            nlsm::NLSM_ID_EPSZ1=parFile["NLSM_ID_EPSZ1"];
            
            nlsm::NLSM_ID_EPSX2=parFile["NLSM_ID_EPSX2"];
            nlsm::NLSM_ID_EPSY2=parFile["NLSM_ID_EPSY2"];
            nlsm::NLSM_ID_EPSZ2=parFile["NLSM_ID_EPSZ2"];

            nlsm::NLSM_ID_R1=parFile["NLSM_ID_R1"];
            nlsm::NLSM_ID_R2=parFile["NLSM_ID_R2"];
            nlsm::NLSM_ID_NU1=parFile["NLSM_ID_NU1"];
            nlsm::NLSM_ID_NU2=parFile["NLSM_ID_NU2"];
            nlsm::NLSM_ID_OMEGA=parFile["NLSM_ID_OMEGA"];

            nlsm::NLSM_WAVELET_TOL=parFile["NLSM_WAVELET_TOL"];
            nlsm::NLSM_CFL_FACTOR=parFile["NLSM_CFL_FACTOR"];

            if(parFile.find("NLSM_MINDEPTH")!=parFile.end())
                nlsm::NLSM_MINDEPTH = parFile["NLSM_MINDEPTH"];

            if(parFile.find("NLSM_WAVE_SPEED_X")!=parFile.end())
                nlsm::NLSM_WAVE_SPEED_X = parFile["NLSM_WAVE_SPEED_X"];
            
            if(parFile.find("NLSM_WAVE_SPEED_Y")!=parFile.end())
                nlsm::NLSM_WAVE_SPEED_Y = parFile["NLSM_WAVE_SPEED_Y"];

            if(parFile.find("NLSM_WAVE_SPEED_Z")!=parFile.end())
                nlsm::NLSM_WAVE_SPEED_Z = parFile["NLSM_WAVE_SPEED_Z"];

            if(parFile.find("NLSM_CHI_REFINE_VAL")!=parFile.end())
                nlsm::NLSM_CHI_REFINE_VAL = parFile["NLSM_CHI_REFINE_VAL"];

            if(parFile.find("NLSM_CHI_COARSEN_VAL")!=parFile.end())
                nlsm::NLSM_CHI_COARSEN_VAL = parFile["NLSM_CHI_COARSEN_VAL"];

            if(parFile.find("NLSM_REFINE_MODE")!=parFile.end())
                nlsm::NLSM_REFINE_MODE = static_cast<nlsm::RefineMode>(parFile["NLSM_REFINE_MODE"]);

            nlsm::NLSM_NUM_REFINE_VARS=parFile["NLSM_NUM_REFINE_VARS"];

            for(unsigned int i=0;i<nlsm::NLSM_NUM_REFINE_VARS;i++)
                nlsm::NLSM_REFINE_VARIABLE_INDICES[i]=parFile["NLSM_REFINE_VARIABLE_INDICES"][i];

            nlsm::NLSM_NUM_EVOL_VARS_VTU_OUTPUT=parFile["NLSM_NUM_EVOL_VARS_VTU_OUTPUT"];

            for(unsigned int i=0;i<nlsm::NLSM_NUM_EVOL_VARS_VTU_OUTPUT;i++)
                nlsm::NLSM_VTU_OUTPUT_EVOL_INDICES[i]=parFile["NLSM_VTU_OUTPUT_EVOL_INDICES"][i];


            vtu_len=NLSM_VTU_FILE_PREFIX.size();
            chp_len=NLSM_CHKPT_FILE_PREFIX.size();
            prf_len=NLSM_PROFILE_FILE_PREFIX.size();

        }

        par::Mpi_Bcast(&NLSM_ELE_ORDER,1,0,comm);
        par::Mpi_Bcast(&NLSM_IO_OUTPUT_FREQ,1,0,comm);
        par::Mpi_Bcast(&NLSM_REMESH_TEST_FREQ,1,0,comm);
        par::Mpi_Bcast(&NLSM_CHECKPT_FREQ,1,0,comm);
        par::Mpi_Bcast(&NLSM_IO_OUTPUT_GAP,1,0,comm);

        par::Mpi_Bcast(&vtu_len,1,0,comm);
        par::Mpi_Bcast(&chp_len,1,0,comm);
        par::Mpi_Bcast(&prf_len,1,0,comm);

        par::Mpi_Bcast(&NLSM_DENDRO_GRAIN_SZ,1,0,comm);
        par::Mpi_Bcast(&NLSM_DENDRO_AMR_FAC,1,0,comm);
        par::Mpi_Bcast(&NLSM_ASYNC_COMM_K,1,0,comm);

        par::Mpi_Bcast(&NLSM_CFL_FACTOR,1,0,comm);

        par::Mpi_Bcast(&NLSM_WAVE_SPEED_X,1,0,comm);
        par::Mpi_Bcast(&NLSM_WAVE_SPEED_Y,1,0,comm);
        par::Mpi_Bcast(&NLSM_WAVE_SPEED_Z,1,0,comm);

        par::Mpi_Bcast(&NLSM_CHI_REFINE_VAL,1,0,comm);
        par::Mpi_Bcast(&NLSM_CHI_COARSEN_VAL,1,0,comm);
        par::Mpi_Bcast((unsigned int*)&NLSM_REFINE_MODE,1,0,comm);
        

        char vtu_name[vtu_len+1];
        char chp_name[chp_len+1];
        char prf_name[prf_len+1];


        if(!rank)
        {
           for(unsigned int k=0;k<vtu_len;k++)
               vtu_name[k]=NLSM_VTU_FILE_PREFIX[k];

            for(unsigned int k=0;k<chp_len;k++)
                chp_name[k]=NLSM_CHKPT_FILE_PREFIX[k];

            for(unsigned int k=0;k<prf_len;k++)
                prf_name[k]=NLSM_PROFILE_FILE_PREFIX[k];

            vtu_name[vtu_len]='\0';
            chp_name[chp_len]='\0';
            prf_name[prf_len]='\0';

        }


        MPI_Bcast(vtu_name,vtu_len+1,MPI_CHAR,0,comm);
        MPI_Bcast(chp_name,chp_len+1,MPI_CHAR,0,comm);
        MPI_Bcast(prf_name,prf_len+1,MPI_CHAR,0,comm);

        NLSM_VTU_FILE_PREFIX=std::string(vtu_name);
        NLSM_CHKPT_FILE_PREFIX=std::string(chp_name);
        NLSM_PROFILE_FILE_PREFIX=std::string(prf_name);


        par::Mpi_Bcast(&NLSM_RESTORE_SOLVER,1,0,comm);
        par::Mpi_Bcast(&NLSM_ENABLE_BLOCK_ADAPTIVITY,1,0,comm);

        par::Mpi_Bcast(&NLSM_WAVELET_TOL,1,0,comm);

        par::Mpi_Bcast(&NLSM_LOAD_IMB_TOL,1,0,comm);
        par::Mpi_Bcast(&NLSM_RK45_TIME_BEGIN,1,0,comm);
        par::Mpi_Bcast(&NLSM_RK45_TIME_END,1,0,comm);
        par::Mpi_Bcast(&NLSM_RK45_TIME_STEP_SIZE,1,0,comm);
        par::Mpi_Bcast(&NLSM_RK45_DESIRED_TOL,1,0,comm);
        par::Mpi_Bcast(&NLSM_DIM,1,0,comm);
        par::Mpi_Bcast(&NLSM_MAXDEPTH,1,0,comm);
        par::Mpi_Bcast(&NLSM_MINDEPTH,1,0,comm);

        par::Mpi_Bcast(&NLSM_ID_TYPE,1,0,comm);
        par::Mpi_Bcast(&NLSM_ID_AMP1,1,0,comm);
        par::Mpi_Bcast(&NLSM_ID_AMP2,1,0,comm);
        par::Mpi_Bcast(&NLSM_ID_DELTA1,1,0,comm);
        par::Mpi_Bcast(&NLSM_ID_DELTA2,1,0,comm);
        par::Mpi_Bcast(&NLSM_ID_XC1,1,0,comm);
        par::Mpi_Bcast(&NLSM_ID_YC1,1,0,comm);
        par::Mpi_Bcast(&NLSM_ID_ZC1,1,0,comm);
        par::Mpi_Bcast(&NLSM_ID_XC2,1,0,comm);
        par::Mpi_Bcast(&NLSM_ID_YC2,1,0,comm);
        par::Mpi_Bcast(&NLSM_ID_ZC2,1,0,comm);
        
        par::Mpi_Bcast(&NLSM_ID_EPSX1,1,0,comm);
        par::Mpi_Bcast(&NLSM_ID_EPSY1,1,0,comm);
        par::Mpi_Bcast(&NLSM_ID_EPSZ1,1,0,comm);

        par::Mpi_Bcast(&NLSM_ID_EPSX2,1,0,comm);
        par::Mpi_Bcast(&NLSM_ID_EPSY2,1,0,comm);
        par::Mpi_Bcast(&NLSM_ID_EPSZ2,1,0,comm);

        par::Mpi_Bcast(&NLSM_ID_R1,1,0,comm);
        par::Mpi_Bcast(&NLSM_ID_R2,1,0,comm);
        par::Mpi_Bcast(&NLSM_ID_NU1,1,0,comm);
        par::Mpi_Bcast(&NLSM_ID_NU2,1,0,comm);
        par::Mpi_Bcast(&NLSM_ID_OMEGA,1,0,comm);

        par::Mpi_Bcast(&NLSM_GRID_MIN_X,1,0,comm);
        par::Mpi_Bcast(&NLSM_GRID_MAX_X,1,0,comm);
        par::Mpi_Bcast(&NLSM_GRID_MIN_Y,1,0,comm);
        par::Mpi_Bcast(&NLSM_GRID_MAX_Y,1,0,comm);
        par::Mpi_Bcast(&NLSM_GRID_MIN_Z,1,0,comm);
        par::Mpi_Bcast(&NLSM_GRID_MAX_Z,1,0,comm);

        par::Mpi_Bcast(&NLSM_BLK_MIN_X,1,0,comm);
        par::Mpi_Bcast(&NLSM_BLK_MIN_Y,1,0,comm);
        par::Mpi_Bcast(&NLSM_BLK_MIN_Z,1,0,comm);

        par::Mpi_Bcast(&NLSM_BLK_MAX_X,1,0,comm);
        par::Mpi_Bcast(&NLSM_BLK_MAX_Y,1,0,comm);
        par::Mpi_Bcast(&NLSM_BLK_MAX_Z,1,0,comm);

        NLSM_OCTREE_MAX[0]=(double )(1u<<nlsm::NLSM_MAXDEPTH);
        NLSM_OCTREE_MAX[1]=(double )(1u<<nlsm::NLSM_MAXDEPTH);
        NLSM_OCTREE_MAX[2]=(double )(1u<<nlsm::NLSM_MAXDEPTH);

        NLSM_COMPD_MIN[0]=NLSM_GRID_MIN_X;
        NLSM_COMPD_MIN[1]=NLSM_GRID_MIN_Y;
        NLSM_COMPD_MIN[2]=NLSM_GRID_MIN_Z;

        NLSM_COMPD_MAX[0]=NLSM_GRID_MAX_X;
        NLSM_COMPD_MAX[1]=NLSM_GRID_MAX_Y;
        NLSM_COMPD_MAX[2]=NLSM_GRID_MAX_Z;

        NLSM_PADDING_WIDTH = NLSM_ELE_ORDER>>1u;
        
        par::Mpi_Bcast(&KO_DISS_SIGMA, 1, 0, comm);

        par::Mpi_Bcast(&NLSM_NUM_REFINE_VARS,1,0,comm);
        par::Mpi_Bcast(&NLSM_NUM_EVOL_VARS_VTU_OUTPUT,1,0,comm);


        if(NLSM_NUM_REFINE_VARS>NLSM_NUM_VARS){std::cout<<"Error[parameter file]: Number of refine variables should be less than number of NLSM_NUM_VARS"<<std::endl; exit(0);}
        if(NLSM_NUM_EVOL_VARS_VTU_OUTPUT>NLSM_NUM_VARS){std::cout<<"Error[parameter file]: Number of evolution VTU variables should be less than number of NLSM_NUM_VARS"<<std::endl; exit(0);}

        par::Mpi_Bcast(NLSM_REFINE_VARIABLE_INDICES,NLSM_NUM_VARS,0,comm);
        par::Mpi_Bcast(NLSM_VTU_OUTPUT_EVOL_INDICES,NLSM_NUM_VARS,0,comm);
        

        

    }

    void dumpParamFile(std::ostream& sout, int root, MPI_Comm comm)
    {
        int rank,npes;
        MPI_Comm_rank(comm,&rank);
        MPI_Comm_size(comm,&npes);

        if(rank==root|| npes==1)
        {
            sout<<"parameters read: "<<std::endl;
            sout<<YLW<<"\tnpes :"<<npes<<NRM<<std::endl;
            sout<<YLW<<"\tNLSM_ELE_ORDER :"<<nlsm::NLSM_ELE_ORDER<<NRM<<std::endl;
            sout<<YLW<<"\tNLSM_PADDING_WIDTH :"<<nlsm::NLSM_PADDING_WIDTH<<NRM<<std::endl;
            sout<<YLW<<"\tNLSM_DIM :"<<nlsm::NLSM_DIM<<NRM<<std::endl;
            sout<<YLW<<"\tNLSM_IO_OUTPUT_FREQ :"<<nlsm::NLSM_IO_OUTPUT_FREQ<<NRM<<std::endl;
            sout<<YLW<<"\tNLSM_REMESH_TEST_FREQ :"<<nlsm::NLSM_REMESH_TEST_FREQ<<NRM<<std::endl;
            sout<<YLW<<"\tNLSM_CHECKPT_FREQ :"<<nlsm::NLSM_CHECKPT_FREQ<<NRM<<std::endl;
            sout<<YLW<<"\tNLSM_RESTORE_SOLVER :"<<nlsm::NLSM_RESTORE_SOLVER<<NRM<<std::endl;
            sout<<YLW<<"\tNLSM_ENABLE_BLOCK_ADAPTIVITY :"<<nlsm::NLSM_ENABLE_BLOCK_ADAPTIVITY<<NRM<<std::endl;
            sout<<YLW<<"\tNLSM_VTU_FILE_PREFIX :"<<nlsm::NLSM_VTU_FILE_PREFIX<<NRM<<std::endl;
            sout<<YLW<<"\tNLSM_CHKPT_FILE_PREFIX :"<<nlsm::NLSM_CHKPT_FILE_PREFIX<<NRM<<std::endl;
            sout<<YLW<<"\tNLSM_PROFILE_FILE_PREFIX :"<<nlsm::NLSM_PROFILE_FILE_PREFIX<<NRM<<std::endl;
            sout<<YLW<<"\tNLSM_IO_OUTPUT_GAP :"<<nlsm::NLSM_IO_OUTPUT_GAP<<NRM<<std::endl;
            sout<<YLW<<"\tNLSM_DENDRO_GRAIN_SZ :"<<nlsm::NLSM_DENDRO_GRAIN_SZ<<NRM<<std::endl;
            sout<<YLW<<"\tNLSM_ASYNC_COMM_K :"<<nlsm::NLSM_ASYNC_COMM_K<<NRM<<std::endl;
            sout<<YLW<<"\tNLSM_DENDRO_AMR_FAC :"<<nlsm::NLSM_DENDRO_AMR_FAC<<NRM<<std::endl;
            sout<<YLW<<"\tNLSM_CFL_FACTOR:"<<nlsm::NLSM_CFL_FACTOR<<NRM<<std::endl;
            sout<<YLW<<"\tNLSM_WAVELET_TOL :"<<nlsm::NLSM_WAVELET_TOL<<NRM<<std::endl;
            sout<<YLW<<"\tNLSM_LOAD_IMB_TOL :"<<nlsm::NLSM_LOAD_IMB_TOL<<NRM<<std::endl;
            sout<<YLW<<"\tNLSM_RK45_TIME_BEGIN :"<<nlsm::NLSM_RK45_TIME_BEGIN<<NRM<<std::endl;
            sout<<YLW<<"\tNLSM_RK45_TIME_END :"<<nlsm::NLSM_RK45_TIME_END<<NRM<<std::endl;
            sout<<YLW<<"\tNLSM_RK45_TIME_STEP_SIZE :"<<nlsm::NLSM_RK45_TIME_STEP_SIZE<<NRM<<std::endl;
            sout<<YLW<<"\tNLSM_RK45_DESIRED_TOL :"<<nlsm::NLSM_RK45_DESIRED_TOL<<NRM<<std::endl;
            sout<<YLW<<"\tNLSM_COMPD_MIN : ( :"<<nlsm::NLSM_COMPD_MIN[0]<<" ,"<<nlsm::NLSM_COMPD_MIN[1]<<","<<nlsm::NLSM_COMPD_MIN[2]<<" )"<<NRM<<std::endl;
            sout<<YLW<<"\tNLSM_COMPD_MAX : ( :"<<nlsm::NLSM_COMPD_MAX[0]<<" ,"<<nlsm::NLSM_COMPD_MAX[1]<<","<<nlsm::NLSM_COMPD_MAX[2]<<" )"<<NRM<<std::endl;
            sout<<YLW<<"\tNLSM_BLK_MIN : ( :"<<nlsm::NLSM_BLK_MIN_X<<" ,"<<nlsm::NLSM_BLK_MIN_Y<<","<<nlsm::NLSM_BLK_MIN_Z<<" )"<<NRM<<std::endl;
            sout<<YLW<<"\tNLSM_BLK_MAX : ( :"<<nlsm::NLSM_BLK_MAX_X<<" ,"<<nlsm::NLSM_BLK_MAX_Y<<","<<nlsm::NLSM_BLK_MAX_Z<<" )"<<NRM<<std::endl;
            sout<<YLW<<"\tNLSM_OCTREE_MIN : ( :"<<nlsm::NLSM_OCTREE_MIN[0]<<" ,"<<nlsm::NLSM_OCTREE_MIN[1]<<","<<nlsm::NLSM_OCTREE_MIN[2]<<" )"<<NRM<<std::endl;
            sout<<YLW<<"\tNLSM_OCTREE_MAX : ( :"<<nlsm::NLSM_OCTREE_MAX[0]<<" ,"<<nlsm::NLSM_OCTREE_MAX[1]<<","<<nlsm::NLSM_OCTREE_MAX[2]<<" )"<<NRM<<std::endl;
            sout<<YLW<<"\tKO_DISS_SIGMA :"<<nlsm::KO_DISS_SIGMA<<NRM<<std::endl;
            sout<<YLW<<"\tNLSM_ID_TYPE:"<<nlsm::NLSM_ID_TYPE<<NRM<<std::endl;
            sout<<YLW<<"\tNLSM_ID_AMP1:"<<nlsm::NLSM_ID_AMP1<<NRM<<std::endl;
            sout<<YLW<<"\tNLSM_ID_AMP2:"<<nlsm::NLSM_ID_AMP2<<NRM<<std::endl;
            sout<<YLW<<"\tNLSM_ID_DELTA1:"<<nlsm::NLSM_ID_DELTA1<<NRM<<std::endl;
            sout<<YLW<<"\tNLSM_ID_DELTA2:"<<nlsm::NLSM_ID_DELTA2<<NRM<<std::endl;
            sout<<YLW<<"\tNLSM_ID_XC1:"<<nlsm::NLSM_ID_XC1<<NRM<<std::endl;
            sout<<YLW<<"\tNLSM_ID_YC1:"<<nlsm::NLSM_ID_YC1<<NRM<<std::endl;
            sout<<YLW<<"\tNLSM_ID_ZC1:"<<nlsm::NLSM_ID_ZC1<<NRM<<std::endl;
            sout<<YLW<<"\tNLSM_ID_XC2:"<<nlsm::NLSM_ID_XC2<<NRM<<std::endl;
            sout<<YLW<<"\tNLSM_ID_YC2:"<<nlsm::NLSM_ID_YC2<<NRM<<std::endl;
            sout<<YLW<<"\tNLSM_ID_ZC2:"<<nlsm::NLSM_ID_ZC2<<NRM<<std::endl;
            sout<<YLW<<"\tNLSM_ID_EPSX1:"<<nlsm::NLSM_ID_EPSX1<<NRM<<std::endl;
            sout<<YLW<<"\tNLSM_ID_EPSY1:"<<nlsm::NLSM_ID_EPSY1<<NRM<<std::endl;
            sout<<YLW<<"\tNLSM_ID_EPSZ1:"<<nlsm::NLSM_ID_EPSY1<<NRM<<std::endl;
            sout<<YLW<<"\tNLSM_ID_EPSX2:"<<nlsm::NLSM_ID_EPSX2<<NRM<<std::endl;
            sout<<YLW<<"\tNLSM_ID_EPSY2:"<<nlsm::NLSM_ID_EPSY2<<NRM<<std::endl;
            sout<<YLW<<"\tNLSM_ID_EPSZ2:"<<nlsm::NLSM_ID_EPSY2<<NRM<<std::endl;
            sout<<YLW<<"\tNLSM_ID_R1:"<<nlsm::NLSM_ID_R1<<NRM<<std::endl;
            sout<<YLW<<"\tNLSM_ID_R2:"<<nlsm::NLSM_ID_R2<<NRM<<std::endl;
            sout<<YLW<<"\tNLSM_ID_NU1:"<<nlsm::NLSM_ID_NU1<<NRM<<std::endl;
            sout<<YLW<<"\tNLSM_ID_NU2:"<<nlsm::NLSM_ID_NU2<<NRM<<std::endl;
            sout<<YLW<<"\tNLSM_ID_OMEGA:"<<nlsm::NLSM_ID_OMEGA<<NRM<<std::endl;
            
            sout<<YLW<<"\tNLSM_DIM :"<<nlsm::NLSM_DIM<<NRM<<std::endl;
            sout<<YLW<<"\tNLSM_MAXDEPTH :"<<nlsm::NLSM_MAXDEPTH<<NRM<<std::endl;
            sout<<YLW<<"\tNLSM_MINDEPTH :"<<nlsm::NLSM_MINDEPTH<<NRM<<std::endl;

            sout<<YLW<<"\tNLSM_NUM_REFINE_VARS :"<<nlsm::NLSM_NUM_REFINE_VARS<<NRM<<std::endl;
            sout<<YLW<<"\tNLSM_REFINE_VARIABLE_INDICES :[";
            for(unsigned int i=0;i<nlsm::NLSM_NUM_REFINE_VARS-1;i++)
                sout<<nlsm::NLSM_REFINE_VARIABLE_INDICES[i]<<", ";
            sout<<nlsm::NLSM_REFINE_VARIABLE_INDICES[nlsm::NLSM_NUM_REFINE_VARS-1]<<"]"<<NRM<<std::endl;

            sout<<YLW<<"\tNLSM_NUM_EVOL_VARS_VTU_OUTPUT :"<<nlsm::NLSM_NUM_EVOL_VARS_VTU_OUTPUT<<NRM<<std::endl;
            sout<<YLW<<"\tNLSM_VTU_OUTPUT_EVOL_INDICES :[";
            for(unsigned int i=0;i<nlsm::NLSM_NUM_EVOL_VARS_VTU_OUTPUT-1;i++)
                sout<<nlsm::NLSM_VTU_OUTPUT_EVOL_INDICES[i]<<", ";
            sout<<nlsm::NLSM_VTU_OUTPUT_EVOL_INDICES[nlsm::NLSM_NUM_EVOL_VARS_VTU_OUTPUT-1]<<"]"<<NRM<<std::endl;

        }

    }

    void initData(const double xx1,const double yy1,const double zz1, double *var)
    {

        const double x=GRIDX_TO_X(xx1);
        const double y=GRIDY_TO_Y(yy1);
        const double z=GRIDZ_TO_Z(zz1);

        const double amp1 = nlsm::NLSM_ID_AMP1;
        const double amp2 = nlsm::NLSM_ID_AMP2;
        const double delta1 = nlsm::NLSM_ID_DELTA1;
        const double delta2 = nlsm::NLSM_ID_DELTA2;
        const double xc1 = nlsm::NLSM_ID_XC1;
        const double yc1 = nlsm::NLSM_ID_YC1;
        const double zc1 = nlsm::NLSM_ID_ZC1;
        const double xc2 = nlsm::NLSM_ID_XC2;
        const double yc2 = nlsm::NLSM_ID_YC2;
        const double zc2 = nlsm::NLSM_ID_ZC2;
        const double epsx1 = nlsm::NLSM_ID_EPSX1;
        const double epsy1 = nlsm::NLSM_ID_EPSY1;
        const double epsz1 = nlsm::NLSM_ID_EPSZ1;
        const double epsx2 = nlsm::NLSM_ID_EPSX2;
        const double epsy2 = nlsm::NLSM_ID_EPSY2;
        const double epsz2 = nlsm::NLSM_ID_EPSZ2;
        const double R1 = nlsm::NLSM_ID_R1;
        const double R2 = nlsm::NLSM_ID_R2;
        const double nu1 = nlsm::NLSM_ID_NU1;
        const double nu2 = nlsm::NLSM_ID_NU2;
        const double Omega = nlsm::NLSM_ID_OMEGA;

        double chi, phi;

        //std::cout<<"initData: "<<x<<", "<<y<<", "<<z<<std::endl;

        #ifdef NLSM_NONLINEAR
            /* regularity requires that chi=0 at the origin for all times. */
                double rsq =  x*x + y*y + z*z;
                if (rsq < 1.0e-13) {
                chi = 0.0;
                phi = 0.0;
                return;
                }
        #endif

        /* if we are not at the origin, then specify a particular ID family */
        if (nlsm::NLSM_ID_TYPE == 0) {
         /* this is the original test data */
          const double amp = amp1;
          const double delta = delta1;
          const double xc = xc1;
          const double yc = yc1;
          const double zc = zc1;
          const double epsx = epsx1;
          const double epsy = epsy1;
          const double epsz = epsz1;
          const double R = R1;
	      double rt = sqrt( epsx*(x-xc)*(x-xc) + epsy*(y-yc)*(y-yc) + epsz*(z-zc)*(z-zc) );
          chi = amp * exp(-(rt-R)*(rt-R)/(delta*delta));
          //chi = amp * exp(-( (x-R)*(x-R) + (y-R)*(y-R) + (z-R)*(z-R) )/(delta*delta));
          phi = 0.0; 

        } else if (nlsm::NLSM_ID_TYPE == 1) {

         /* this is family (a) from Liebling, arXiv:gr-qc/0202093 */

          double rt1 = sqrt( epsx1*(x-xc1)*(x-xc1) + epsy1*(y-yc1)*(y-yc1) + (z-zc1)*(z-zc1) );
          rt1 += 1.0e-14;
          chi = amp1 * exp(-(rt1-R1)*(rt1-R1)/(delta1*delta1));
          double dGdr = amp1 * exp( -(rt1-R1)*(rt1-R1)/(delta1*delta1) ) 
                      * ( -2.0*(rt1-R1)/(delta1*delta1) );
          double dGdx = -2.0*chi*(rt1-R1)/(delta1*delta1)*epsx1*(x-xc1)/rt1;
          double dGdy = -2.0*chi*(rt1-R1)/(delta1*delta1)*epsy1*(y-yc1)/rt1;
          phi = nu1*dGdr + Omega*(y*dGdx - x*dGdy); 
           
        } else if (nlsm::NLSM_ID_TYPE == 2) {

         /* this is family (b) from Liebling, arXiv:gr-qc/0202093 */

          double rt1 = sqrt( epsx1*(x-xc1)*(x-xc1) + epsy1*(y-yc1)*(y-yc1) + (z-zc1)*(z-zc1) );
          rt1 += 1.0e-14;
          double chi1 = amp1 * exp(-(rt1-R1)*(rt1-R1)/(delta1*delta1));

          double rt2 = sqrt( epsx2*(x-xc2)*(x-xc2) + epsy2*(y-yc2)*(y-yc2) + (z-zc2)*(z-zc2) );
          rt2 += 1.0e-14;
          double chi2 = amp2 * exp(-(rt2-R2)*(rt2-R2)/(delta2*delta2));
          double dGdx1 = -2.0*chi1*(rt1-R1)/(delta1*delta1)*epsx1*(x-xc1)/rt1;
          double dGdx2 = -2.0*chi2*(rt2-R2)/(delta2*delta2)*epsx2*(x-xc2)/rt2;
          chi = chi1 + chi2;
          phi = nu1*dGdx1 + nu2*dGdx2;

        } else if (nlsm::NLSM_ID_TYPE == 3) {

        // simplistic initial data to make the comparison with analytical solution easier

          const double amp = amp1;
          const double delta = delta1;
          const double xc = xc1;
          const double yc = yc1;
          const double zc = zc1;
          const double epsx = epsx1;
          const double epsy = epsy1;
          const double epsz = epsz1;
          const double R = R1;
	      double rt = epsx*(x-xc)+epsy*(y-yc)+epsz*(z-zc);
          chi = amp * exp(-(rt-R)/(delta*delta));
          phi =0;//-sqrt(3)*amp * exp(-(rt-R)/(delta*delta)); 

        } else if (nlsm::NLSM_ID_TYPE == 4) {

          double rt1 = sqrt( (x-xc1)*(x-xc1) + (y-yc1)*(y-yc1) + (z-zc1)*(z-zc1) );

          if (rt1 < 1.0e-7) {
            chi = 0.0;
            phi = 0.0;
          }
          else {
            chi = amp1 * exp(-(rt1-R1)*(rt1-R1)/(delta1*delta1)) / rt1;
            phi = chi * ( 2.0*(rt1-R1)/(delta1*delta1) );
          }
          //printf("x=%g, y=%g, z=%g, chi=%g phi=%g\n",x,y,z,chi,phi);

        } else {
          chi = 0.0;
          phi = 0.0;
        }

        var[VAR::U_CHI] = chi;
        var[VAR::U_PHI] = phi;

    }


    void analyticalSol(const double xx1,const double yy1,const double zz1,const double t, double *var)
    {
        const double x=GRIDX_TO_X(xx1);
        const double y=GRIDY_TO_Y(yy1);
        const double z=GRIDZ_TO_Z(zz1);

        const double amp1 = nlsm::NLSM_ID_AMP1;
        const double amp2 = nlsm::NLSM_ID_AMP2;
        const double delta1 = nlsm::NLSM_ID_DELTA1;
        const double delta2 = nlsm::NLSM_ID_DELTA2;
        const double xc1 = nlsm::NLSM_ID_XC1;
        const double yc1 = nlsm::NLSM_ID_YC1;
        const double zc1 = nlsm::NLSM_ID_ZC1;
        const double xc2 = nlsm::NLSM_ID_XC2;
        const double yc2 = nlsm::NLSM_ID_YC2;
        const double zc2 = nlsm::NLSM_ID_ZC2;
        const double epsx1 = nlsm::NLSM_ID_EPSX1;
        const double epsy1 = nlsm::NLSM_ID_EPSY1;
        const double epsz1 = nlsm::NLSM_ID_EPSZ1;
        const double epsx2 = nlsm::NLSM_ID_EPSX2;
        const double epsy2 = nlsm::NLSM_ID_EPSY2;
        const double epsz2 = nlsm::NLSM_ID_EPSZ2;
        const double R1 = nlsm::NLSM_ID_R1;
        const double R2 = nlsm::NLSM_ID_R2;
        const double nu1 = nlsm::NLSM_ID_NU1;
        const double nu2 = nlsm::NLSM_ID_NU2;
        const double Omega = nlsm::NLSM_ID_OMEGA;

        double chi=0.0, phi=0.0;

        if (nlsm::NLSM_ID_TYPE == 0) {
         /* this is the original test data */
          const double amp = amp1;
          const double delta = delta1;
          const double xc = xc1;
          const double yc = yc1;
          const double zc = zc1;
          const double epsx = epsx1;
          const double epsy = epsy1;
          const double epsz = epsz1;
          const double R = R1;
	      double rt = sqrt( epsx*(x-xc)*(x-xc) + epsy*(y-yc)*(y-yc) + epsz*(z-zc)*(z-zc) );
          double rt_pt = sqrt( epsx*(x + nlsm::NLSM_WAVE_SPEED_X*t -xc)*(x + nlsm::NLSM_WAVE_SPEED_X*t -xc) + epsy*(y + nlsm::NLSM_WAVE_SPEED_Y*t -yc)*(y + nlsm::NLSM_WAVE_SPEED_Y*t -yc) + epsz*(z + nlsm::NLSM_WAVE_SPEED_Z*t -zc)*(z + nlsm::NLSM_WAVE_SPEED_Z*t -zc) ); 
          double rt_mt = sqrt( epsx*(x - nlsm::NLSM_WAVE_SPEED_X*t -xc)*(x - nlsm::NLSM_WAVE_SPEED_X*t -xc) + epsy*(y - nlsm::NLSM_WAVE_SPEED_Y*t -yc)*(y - nlsm::NLSM_WAVE_SPEED_Y*t -yc) + epsz*(z - nlsm::NLSM_WAVE_SPEED_Z*t -zc)*(z - nlsm::NLSM_WAVE_SPEED_Z*t -zc) );
          chi = 0.5*amp *(exp(-(rt_pt-R)*(rt_pt-R)/(delta*delta)) + exp(-(rt_mt-R)*(rt_mt-R)/(delta*delta))); 
          phi = 0.0; 

        }else if (nlsm::NLSM_ID_TYPE == 3) {
            /* this is the original test data */
            const double amp = 1.0e-2;
            const double delta = 1.0;
            const double xc = 0.0;
            const double yc = 0.0;
            const double zc = 0.0;
            const double epsx = 1.0;
            const double epsy = 1.0;
            const double epsz = 1.0;
            const double R = 0.0;
            double rt = ( epsx*(x-xc)*(x-xc) + epsy*(y-yc)*(y-yc) + epsz*(z-zc)*(z-zc) );
            double rt_plus_t  = (epsx*(x + t - xc)*(x + t - xc) + epsy*(y + t - yc)*(y + t - yc) + epsz*(z + t - zc)*(z + t - zc) );
            double rt_minus_t = (epsx*(x - t - xc)*(x - t - xc) + epsy*(y - t - yc)*(y - t - yc) + epsz*(z - t - zc)*(z - t - zc) );
            chi = amp *exp(-(rt_plus_t-R)/(delta*delta));
            phi = 0.0;
            
        } else if (nlsm::NLSM_ID_TYPE == 4) {
            double rt1 = sqrt( (x-xc1)*(x-xc1) + (y-yc1)*(y-yc1) + (z-zc1)*(z-zc1) );

          if (rt1 < 1.0e-7) {
            chi = 0.0;
            phi = 0.0;
          }
          else {
            chi = amp1 * exp(-(rt1-R1-t)*(rt1-R1-t)/(delta1*delta1)) / rt1;
            phi = chi * ( 2.0*(rt1-R1-t)/(delta1*delta1) );
          }

        }

        var[VAR::U_CHI] = chi;
        var[VAR::U_PHI] = phi;

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
        xRange_e=pt_g_max[0];//((rank+1)*(pt_g_max[0]-pt_g_min[0]))/npes + pt_g_min[0];

        // xRange_b=(rank*(pt_g_max[0]-pt_g_min[0]))/npes + pt_g_min[0];
        // xRange_e=((rank+1)*(pt_g_max[0]-pt_g_min[0]))/npes + pt_g_min[0];

        unsigned int stepSz=1u<<(maxDepth-regLev);

       /* std::cout<<" x min: "<<xRange_b<<" x_max: "<<xRange_e<<std::endl;
        std::cout<<" y min: "<<yRange_b<<" y_max: "<<yRange_e<<std::endl;
        std::cout<<" z min: "<<zRange_b<<" z_max: "<<zRange_e<<std::endl;*/

        // note that other wise this will create lot of duplicates. 
        if(!rank)
        for(unsigned int x=xRange_b;x<xRange_e;x+=stepSz)
            for(unsigned int y=yRange_b;y<yRange_e;y+=stepSz)
              for(unsigned int z=zRange_b;z<zRange_e;z+=stepSz)
              {
                if(x>=(1u<<maxDepth)) x=x-1;	
		        if(y>=(1u<<maxDepth)) y=y-1;	
		        if(z>=(1u<<maxDepth)) z=z-1;	

                tmpNodes.push_back(ot::TreeNode(x,y,z,regLev,m_uiDim,maxDepth));
              }
        
        // now scatter the elements.
        DendroIntL totalNumOcts = tmpNodes.size(), numOcts;
        par::Mpi_Bcast<DendroIntL>(&totalNumOcts, 1, 0, comm);

        const unsigned int size = npes;
        std::vector<ot::TreeNode> nodes_new;
        // TODO do proper load balancing.
        numOcts = totalNumOcts/size + (rank < totalNumOcts%size);
        par::scatterValues<ot::TreeNode>(tmpNodes, nodes_new, numOcts, comm);
        std::swap(tmpNodes, nodes_new);
        nodes_new.clear();


        return ;


    }

    double computeWTol(double x,double y,double z,double* hx)
    {
       return nlsm::NLSM_WAVELET_TOL;
    }


    bool isRemeshForce(const ot::Mesh* pMesh, const double ** unzipVec, unsigned int vIndex, double refine_th, double coarsen_th, bool isOverwrite)
    {
        const unsigned int eleLocalBegin = pMesh->getElementLocalBegin();
        const unsigned int eleLocalEnd = pMesh->getElementLocalEnd();
        bool isOctChange=false;
        bool isOctChange_g =false;
        const unsigned int eOrder = pMesh->getElementOrder();

        if(pMesh->isActive())
        {
            ot::TreeNode * pNodes = (ot::TreeNode*) &(*(pMesh->getAllElements().begin()));
            
            if(isOverwrite)
            for(unsigned int ele = eleLocalBegin; ele< eleLocalEnd; ele++)
                pNodes[ele].setFlag(((OCT_NO_CHANGE<<NUM_LEVEL_BITS)|pNodes[ele].getLevel()));


            const std::vector<ot::Block>& blkList = pMesh->getLocalBlockList();
            unsigned int sz[3];
            unsigned int ei[3];
            
            // refine test
            for(unsigned int b=0; b< blkList.size(); b++)
            {
                const ot::TreeNode blkNode = blkList[b].getBlockNode();

                sz[0]=blkList[b].getAllocationSzX();
                sz[1]=blkList[b].getAllocationSzY();
                sz[2]=blkList[b].getAllocationSzZ();

                const unsigned int bflag = blkList[b].getBlkNodeFlag();
                const unsigned int offset = blkList[b].getOffset();

                const unsigned int regLev=blkList[b].getRegularGridLev();
                const unsigned int eleIndexMax=(1u<<(regLev-blkNode.getLevel()))-1;
                const unsigned int eleIndexMin=0;

                for(unsigned int ele = blkList[b].getLocalElementBegin(); ele< blkList[b].getLocalElementEnd(); ele++)
                {
                    ei[0]=(pNodes[ele].getX()-blkNode.getX())>>(m_uiMaxDepth-regLev);
                    ei[1]=(pNodes[ele].getY()-blkNode.getY())>>(m_uiMaxDepth-regLev);
                    ei[2]=(pNodes[ele].getZ()-blkNode.getZ())>>(m_uiMaxDepth-regLev);

                    if((bflag &(1u<<OCT_DIR_LEFT)) && ei[0]==eleIndexMin)   continue;
                    if((bflag &(1u<<OCT_DIR_DOWN)) && ei[1]==eleIndexMin)   continue;
                    if((bflag &(1u<<OCT_DIR_BACK)) && ei[2]==eleIndexMin)   continue;

                    if((bflag &(1u<<OCT_DIR_RIGHT)) && ei[0]==eleIndexMax)  continue;
                    if((bflag &(1u<<OCT_DIR_UP)) && ei[1]==eleIndexMax)     continue;
                    if((bflag &(1u<<OCT_DIR_FRONT)) && ei[2]==eleIndexMax)  continue;

                    // refine test. 
                    for(unsigned int k=3; k< eOrder+1 +   3; k++)
                     for(unsigned int j=3; j< eOrder+1 +  3; j++)
                      for(unsigned int i=3; i< eOrder+1 + 3; i++)
                      {
                          if ( (unzipVec[vIndex][offset + (ei[2]*eOrder + k)*sz[0]*sz[1] + (ei[1]*eOrder + j)*sz[0] + (ei[0]*eOrder + i)] < refine_th) &&  (unzipVec[vIndex][offset + (ei[2]*eOrder + k)*sz[0]*sz[1] + (ei[1]*eOrder + j)*sz[0] + (ei[0]*eOrder + i)]) > coarsen_th )
                          {
                            if( (pNodes[ele].getLevel() + MAXDEAPTH_LEVEL_DIFF +1) < m_uiMaxDepth  )
                                pNodes[ele].setFlag(((OCT_SPLIT<<NUM_LEVEL_BITS)|pNodes[ele].getLevel()));

                          }

                      }
                    
                    
                }

            }

            //coarsen test. 
            for(unsigned int b=0; b< blkList.size(); b++)
            {
                const ot::TreeNode blkNode = blkList[b].getBlockNode();

                sz[0]=blkList[b].getAllocationSzX();
                sz[1]=blkList[b].getAllocationSzY();
                sz[2]=blkList[b].getAllocationSzZ();

                const unsigned int bflag = blkList[b].getBlkNodeFlag();
                const unsigned int offset = blkList[b].getOffset();

                const unsigned int regLev=blkList[b].getRegularGridLev();
                const unsigned int eleIndexMax=(1u<<(regLev-blkNode.getLevel()))-1;
                const unsigned int eleIndexMin=0;

                if((eleIndexMax==0) || (bflag!=0)) continue; // this implies the blocks with only 1 child and boundary blocks.

                for(unsigned int ele = blkList[b].getLocalElementBegin(); ele< blkList[b].getLocalElementEnd(); ele++)
                {
                    assert(pNodes[ele].getParent()==pNodes[ele+NUM_CHILDREN-1].getParent());
                    bool isCoarsen =true;

                    for(unsigned int child=0;child<NUM_CHILDREN;child++)
                    {
                        if((pNodes[ele+child].getFlag()>>NUM_LEVEL_BITS)==OCT_SPLIT)
                        {
                            isCoarsen=false;
                            break;
                        }

                    }

                    if(isCoarsen && pNodes[ele].getLevel()>1)
                    {
                        bool coarse = true;
                        for(unsigned int child=0;child<NUM_CHILDREN;child++)
                        {
                            ei[0]=(pNodes[ele + child].getX()-blkNode.getX())>>(m_uiMaxDepth-regLev);
                            ei[1]=(pNodes[ele + child].getY()-blkNode.getY())>>(m_uiMaxDepth-regLev);
                            ei[2]=(pNodes[ele + child].getZ()-blkNode.getZ())>>(m_uiMaxDepth-regLev);

                            for(unsigned int k=3; k< eOrder+1 + 3; k++)
                            for(unsigned int j=3; j< eOrder+1 +3; j++)
                             for(unsigned int i=3; i< eOrder+ +3; i++)
                             {
                                if ( (unzipVec[vIndex][offset + (ei[2]*eOrder + k)*sz[0]*sz[1] + (ei[1]*eOrder + j)*sz[0] + (ei[0]*eOrder + i)]) > coarsen_th )
                                    coarse = false;
                             }
                        }

                        if(coarse)
                            for(unsigned int child=0;child<NUM_CHILDREN;child++)
                                pNodes[ele+child].setFlag(((OCT_COARSE<<NUM_LEVEL_BITS)|pNodes[ele].getLevel()));


                    }

                    ele = ele + NUM_CHILDREN-1;

                    
                    
                    
                }

            }



            for(unsigned int ele=eleLocalBegin;ele<eleLocalEnd;ele++)
                if((pNodes[ele].getFlag()>>NUM_LEVEL_BITS)==OCT_SPLIT) // trigger remesh only when some refinement occurs (laid back remesh :)  )
                { 
                    isOctChange=true;
                    break;
                }

            

        }

        bool isOctChanged_g;
        MPI_Allreduce(&isOctChange,&isOctChanged_g,1,MPI_CXX_BOOL,MPI_LOR,pMesh->getMPIGlobalCommunicator());
        //if(!m_uiGlobalRank) std::cout<<"is oct changed: "<<isOctChanged_g<<std::endl;
        return isOctChanged_g;
    }

    unsigned int getEleWeight(const ot::TreeNode* pNode) 
    {
        return (1u<<(3*pNode->getLevel()))*1;
    }

}// end of namespace nlsm





namespace nlsm
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

            t_rhs_a.start();
            t_rhs_b.start();
            t_rhs_gt.start();
            t_rhs_chi.start();
            t_rhs_At.start();
            t_rhs_K.start();
            t_rhs_Gt.start();
            t_rhs_B.start();


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

            t_rhs_a.snapreset();
            t_rhs_b.snapreset();
            t_rhs_gt.snapreset();
            t_rhs_chi.snapreset();
            t_rhs_At.snapreset();
            t_rhs_K.snapreset();
            t_rhs_Gt.snapreset();
            t_rhs_B.snapreset();

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
                outfile<<"partition tol : "<<nlsm::NLSM_LOAD_IMB_TOL<<std::endl;
                outfile<<"wavelet tol : "<<nlsm::NLSM_WAVELET_TOL<<std::endl;
                outfile<<"maxdepth : "<<nlsm::NLSM_MAXDEPTH<<std::endl;

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


#ifdef NLSM_PROFILE_HUMAN_READABLE
            if(!activeRank)
            {
                sprintf(fName,"%s_im.prof",filePrefix);
                outfile.open (fName,std::fstream::app);
                if(outfile.fail()) {std::cout<<fName<<" file open failed "<<std::endl; return ;}

                outfile<<"active npes : "<<activeNpes<<std::endl;
                outfile<<"global npes : "<<globalNpes<<std::endl;
                outfile<<"current step : "<<currentStep<<std::endl;
                outfile<<"partition tol : "<<nlsm::NLSM_LOAD_IMB_TOL<<std::endl;
                outfile<<"wavelet tol : "<<nlsm::NLSM_WAVELET_TOL<<std::endl;
                outfile<<"maxdepth : "<<nlsm::NLSM_MAXDEPTH<<std::endl;

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



            t_stat=t_rhs_a.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --compute_rhs_a ";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=t_rhs_b.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --compute_rhs_b ";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=t_rhs_gt.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --compute_rhs_gt ";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_rhs_chi.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --compute_rhs_chi ";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=t_rhs_At.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --compute_rhs_At ";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=t_rhs_K.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --compute_rhs_K ";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=t_rhs_Gt.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --compute_rhs_Gt ";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=t_rhs_B.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --compute_rhs_B ";
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
            if(!rank) outfile<<nlsm::NLSM_LOAD_IMB_TOL<<"\t ";
            if(!rank) outfile<<nlsm::NLSM_WAVELET_TOL<<"\t ";
            if(!rank) outfile<<nlsm::NLSM_MAXDEPTH<<"\t ";

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
