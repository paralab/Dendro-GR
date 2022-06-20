//
// Created by milinda on 7/26/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief Contains utility functions for EM3 simulation.
*/
//

#include "em3Utils.h"

namespace em3
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

            if(parFile.find("EM3_ELE_ORDER")!=parFile.end())
                em3::EM3_ELE_ORDER = parFile["EM3_ELE_ORDER"];

            em3::EM3_IO_OUTPUT_FREQ=parFile["EM3_IO_OUTPUT_FREQ"];
            em3::EM3_REMESH_TEST_FREQ=parFile["EM3_REMESH_TEST_FREQ"];
            em3::EM3_CHECKPT_FREQ=parFile["EM3_CHECKPT_FREQ"];
            em3::EM3_IO_OUTPUT_GAP=parFile["EM3_IO_OUTPUT_GAP"];
            em3::EM3_VTU_FILE_PREFIX=parFile["EM3_VTU_FILE_PREFIX"].get<std::string>();
            em3::EM3_CHKPT_FILE_PREFIX=parFile["EM3_CHKPT_FILE_PREFIX"].get<std::string>();
            em3::EM3_PROFILE_FILE_PREFIX=parFile["EM3_PROFILE_FILE_PREFIX"].get<std::string>();
            em3::EM3_RESTORE_SOLVER=parFile["EM3_RESTORE_SOLVER"];
            em3::EM3_ENABLE_BLOCK_ADAPTIVITY=parFile["EM3_ENABLE_BLOCK_ADAPTIVITY"];
            em3::EM3_BLK_MIN_X=parFile["EM3_BLK_MIN_X"];
            em3::EM3_BLK_MIN_Y=parFile["EM3_BLK_MIN_Y"];
            em3::EM3_BLK_MIN_Z=parFile["EM3_BLK_MIN_Z"];
            em3::EM3_BLK_MAX_X=parFile["EM3_BLK_MAX_X"];
            em3::EM3_BLK_MAX_Y=parFile["EM3_BLK_MAX_Y"];
            em3::EM3_BLK_MAX_Z=parFile["EM3_BLK_MAX_Z"];
            em3::EM3_DENDRO_GRAIN_SZ=parFile["EM3_DENDRO_GRAIN_SZ"];
            em3::EM3_ASYNC_COMM_K=parFile["EM3_ASYNC_COMM_K"];
            em3::EM3_DENDRO_AMR_FAC=parFile["EM3_DENDRO_AMR_FAC"];
            em3::EM3_LOAD_IMB_TOL=parFile["EM3_LOAD_IMB_TOL"];
            em3::EM3_RK45_TIME_BEGIN=parFile["EM3_RK45_TIME_BEGIN"];
            em3::EM3_RK45_TIME_END=parFile["EM3_RK45_TIME_END"];
            em3::EM3_RK45_TIME_STEP_SIZE=parFile["EM3_RK45_TIME_STEP_SIZE"];
            em3::EM3_RK45_DESIRED_TOL=parFile["EM3_RK45_DESIRED_TOL"];
            em3::EM3_DIM=parFile["EM3_DIM"];
            em3::EM3_MAXDEPTH=parFile["EM3_MAXDEPTH"];
            em3::EM3_GRID_MIN_X=parFile["EM3_GRID_MIN_X"];
            em3::EM3_GRID_MAX_X=parFile["EM3_GRID_MAX_X"];
            em3::EM3_GRID_MIN_Y=parFile["EM3_GRID_MIN_Y"];
            em3::EM3_GRID_MAX_Y=parFile["EM3_GRID_MAX_Y"];
            em3::EM3_GRID_MIN_Z=parFile["EM3_GRID_MIN_Z"];
            em3::EM3_GRID_MAX_Z=parFile["EM3_GRID_MAX_Z"];
            em3::KO_DISS_SIGMA=parFile["KO_DISS_SIGMA"];
            if(parFile.find("DISSIPATION_TYPE") != parFile.end()) {
                em3::DISSIPATION_TYPE=parFile["DISSIPATION_TYPE"];
            };

            em3::EM3_ID_TYPE=parFile["EM3_ID_TYPE"];
            em3::EM3_ID_AMP1=parFile["EM3_ID_AMP1"];
            em3::EM3_ID_LAMBDA1=parFile["EM3_ID_LAMBDA1"];
            em3::EM3_ID_AMP2=parFile["EM3_ID_AMP2"];
            //em3::EM3_ID_DELTA1=parFile["EM3_ID_DELTA1"];
            //em3::EM3_ID_DELTA2=parFile["EM3_ID_DELTA2"];
            //em3::EM3_ID_XC1=parFile["EM3_ID_XC1"];
            //em3::EM3_ID_YC1=parFile["EM3_ID_YC1"];
            //em3::EM3_ID_ZC1=parFile["EM3_ID_ZC1"];
            //em3::EM3_ID_XC2=parFile["EM3_ID_XC2"];
            //em3::EM3_ID_YC2=parFile["EM3_ID_YC2"];
            //em3::EM3_ID_ZC2=parFile["EM3_ID_ZC2"];
            
            //em3::EM3_ID_EPSX1=parFile["EM3_ID_EPSX1"];
            //em3::EM3_ID_EPSY1=parFile["EM3_ID_EPSY1"];
            //em3::EM3_ID_EPSZ1=parFile["EM3_ID_EPSZ1"];
            
            //em3::EM3_ID_EPSX2=parFile["EM3_ID_EPSX2"];
            //em3::EM3_ID_EPSY2=parFile["EM3_ID_EPSY2"];
            //em3::EM3_ID_EPSZ2=parFile["EM3_ID_EPSZ2"];

            //em3::EM3_ID_R1=parFile["EM3_ID_R1"];
            //em3::EM3_ID_R2=parFile["EM3_ID_R2"];
            //em3::EM3_ID_NU1=parFile["EM3_ID_NU1"];
            //em3::EM3_ID_NU2=parFile["EM3_ID_NU2"];
            //em3::EM3_ID_OMEGA=parFile["EM3_ID_OMEGA"];

            em3::EM3_WAVELET_TOL=parFile["EM3_WAVELET_TOL"];
            em3::EM3_CFL_FACTOR=parFile["EM3_CFL_FACTOR"];


            if(parFile.find("EM3_CHI_REFINE_VAL")!=parFile.end())
                em3::EM3_CHI_REFINE_VAL = parFile["EM3_CHI_REFINE_VAL"];

            if(parFile.find("EM3_VTU_Z_SLICE_ONLY")!=parFile.end())
                em3::EM3_VTU_Z_SLICE_ONLY = parFile["EM3_VTU_Z_SLICE_ONLY"];

            if(parFile.find("EM3_CHI_COARSEN_VAL")!=parFile.end())
                em3::EM3_CHI_COARSEN_VAL = parFile["EM3_CHI_COARSEN_VAL"];

            if(parFile.find("EM3_REFINEMENT_MODE")!=parFile.end())
                em3::EM3_REFINEMENT_MODE = static_cast<em3::RefinementMode>(parFile["EM3_REFINEMENT_MODE"]);

            em3::EM3_NUM_REFINE_VARS=parFile["EM3_NUM_REFINE_VARS"];

            for(unsigned int i=0;i<em3::EM3_NUM_REFINE_VARS;i++)
                em3::EM3_REFINE_VARIABLE_INDICES[i]=parFile["EM3_REFINE_VARIABLE_INDICES"][i];

            em3::EM3_NUM_EVOL_VARS_VTU_OUTPUT=parFile["EM3_NUM_EVOL_VARS_VTU_OUTPUT"];
            em3::EM3_NUM_CONST_VARS_VTU_OUTPUT=parFile["EM3_NUM_CONST_VARS_VTU_OUTPUT"];

            for(unsigned int i=0;i<em3::EM3_NUM_EVOL_VARS_VTU_OUTPUT;i++)
                em3::EM3_VTU_OUTPUT_EVOL_INDICES[i]=parFile["EM3_VTU_OUTPUT_EVOL_INDICES"][i];

            for(unsigned int i=0;i<em3::EM3_NUM_CONST_VARS_VTU_OUTPUT;i++)
                em3::EM3_VTU_OUTPUT_CONST_INDICES[i]=parFile["EM3_VTU_OUTPUT_CONST_INDICES"][i];

            vtu_len=EM3_VTU_FILE_PREFIX.size();
            chp_len=EM3_CHKPT_FILE_PREFIX.size();
            prf_len=EM3_PROFILE_FILE_PREFIX.size();

        }

        par::Mpi_Bcast(&EM3_ELE_ORDER,1,0,comm);
        par::Mpi_Bcast(&EM3_IO_OUTPUT_FREQ,1,0,comm);
        par::Mpi_Bcast(&EM3_REMESH_TEST_FREQ,1,0,comm);
        par::Mpi_Bcast(&EM3_CHECKPT_FREQ,1,0,comm);
        par::Mpi_Bcast(&EM3_IO_OUTPUT_GAP,1,0,comm);

        par::Mpi_Bcast(&vtu_len,1,0,comm);
        par::Mpi_Bcast(&chp_len,1,0,comm);
        par::Mpi_Bcast(&prf_len,1,0,comm);

        par::Mpi_Bcast(&EM3_DENDRO_GRAIN_SZ,1,0,comm);
        par::Mpi_Bcast(&EM3_DENDRO_AMR_FAC,1,0,comm);
        par::Mpi_Bcast(&EM3_ASYNC_COMM_K,1,0,comm);

        par::Mpi_Bcast(&EM3_CFL_FACTOR,1,0,comm);

        par::Mpi_Bcast(&EM3_CHI_REFINE_VAL,1,0,comm);
        par::Mpi_Bcast(&EM3_CHI_COARSEN_VAL,1,0,comm);
        par::Mpi_Bcast((unsigned int*)&EM3_REFINEMENT_MODE,1,0,comm);
        
        par::Mpi_Bcast(&EM3_VTU_Z_SLICE_ONLY,1,0,comm);


        char vtu_name[vtu_len+1];
        char chp_name[chp_len+1];
        char prf_name[prf_len+1];


        if(!rank)
        {
           for(unsigned int k=0;k<vtu_len;k++)
               vtu_name[k]=EM3_VTU_FILE_PREFIX[k];

            for(unsigned int k=0;k<chp_len;k++)
                chp_name[k]=EM3_CHKPT_FILE_PREFIX[k];

            for(unsigned int k=0;k<prf_len;k++)
                prf_name[k]=EM3_PROFILE_FILE_PREFIX[k];

            vtu_name[vtu_len]='\0';
            chp_name[chp_len]='\0';
            prf_name[prf_len]='\0';

        }


        MPI_Bcast(vtu_name,vtu_len+1,MPI_CHAR,0,comm);
        MPI_Bcast(chp_name,chp_len+1,MPI_CHAR,0,comm);
        MPI_Bcast(prf_name,prf_len+1,MPI_CHAR,0,comm);

        EM3_VTU_FILE_PREFIX=std::string(vtu_name);
        EM3_CHKPT_FILE_PREFIX=std::string(chp_name);
        EM3_PROFILE_FILE_PREFIX=std::string(prf_name);


        par::Mpi_Bcast(&EM3_RESTORE_SOLVER,1,0,comm);
        par::Mpi_Bcast(&EM3_ENABLE_BLOCK_ADAPTIVITY,1,0,comm);

        par::Mpi_Bcast(&EM3_WAVELET_TOL,1,0,comm);

        par::Mpi_Bcast(&EM3_LOAD_IMB_TOL,1,0,comm);
        par::Mpi_Bcast(&EM3_RK45_TIME_BEGIN,1,0,comm);
        par::Mpi_Bcast(&EM3_RK45_TIME_END,1,0,comm);
        par::Mpi_Bcast(&EM3_RK45_TIME_STEP_SIZE,1,0,comm);
        par::Mpi_Bcast(&EM3_RK45_DESIRED_TOL,1,0,comm);
        par::Mpi_Bcast(&EM3_DIM,1,0,comm);
        par::Mpi_Bcast(&EM3_MAXDEPTH,1,0,comm);

        par::Mpi_Bcast(&EM3_ID_TYPE,1,0,comm);
        par::Mpi_Bcast(&EM3_ID_AMP1,1,0,comm);
        par::Mpi_Bcast(&EM3_ID_LAMBDA1,1,0,comm);
        par::Mpi_Bcast(&EM3_ID_AMP2,1,0,comm);
        //par::Mpi_Bcast(&EM3_ID_DELTA1,1,0,comm);
        //par::Mpi_Bcast(&EM3_ID_DELTA2,1,0,comm);
        //par::Mpi_Bcast(&EM3_ID_XC1,1,0,comm);
        //par::Mpi_Bcast(&EM3_ID_YC1,1,0,comm);
        //par::Mpi_Bcast(&EM3_ID_ZC1,1,0,comm);
        //par::Mpi_Bcast(&EM3_ID_XC2,1,0,comm);
        //par::Mpi_Bcast(&EM3_ID_YC2,1,0,comm);
        //par::Mpi_Bcast(&EM3_ID_ZC2,1,0,comm);

        //par::Mpi_Bcast(&EM3_ID_EPSX1,1,0,comm);
        //par::Mpi_Bcast(&EM3_ID_EPSY1,1,0,comm);
        //par::Mpi_Bcast(&EM3_ID_EPSZ1,1,0,comm);

        //par::Mpi_Bcast(&EM3_ID_EPSX2,1,0,comm);
        //par::Mpi_Bcast(&EM3_ID_EPSY2,1,0,comm);
        //par::Mpi_Bcast(&EM3_ID_EPSZ2,1,0,comm);

        //par::Mpi_Bcast(&EM3_ID_R1,1,0,comm);
        //par::Mpi_Bcast(&EM3_ID_R2,1,0,comm);
        //par::Mpi_Bcast(&EM3_ID_NU1,1,0,comm);
        //par::Mpi_Bcast(&EM3_ID_NU2,1,0,comm);
        //par::Mpi_Bcast(&EM3_ID_OMEGA,1,0,comm);

        par::Mpi_Bcast(&EM3_GRID_MIN_X,1,0,comm);
        par::Mpi_Bcast(&EM3_GRID_MAX_X,1,0,comm);
        par::Mpi_Bcast(&EM3_GRID_MIN_Y,1,0,comm);
        par::Mpi_Bcast(&EM3_GRID_MAX_Y,1,0,comm);
        par::Mpi_Bcast(&EM3_GRID_MIN_Z,1,0,comm);
        par::Mpi_Bcast(&EM3_GRID_MAX_Z,1,0,comm);

        par::Mpi_Bcast(&EM3_BLK_MIN_X,1,0,comm);
        par::Mpi_Bcast(&EM3_BLK_MIN_Y,1,0,comm);
        par::Mpi_Bcast(&EM3_BLK_MIN_Z,1,0,comm);

        par::Mpi_Bcast(&EM3_BLK_MAX_X,1,0,comm);
        par::Mpi_Bcast(&EM3_BLK_MAX_Y,1,0,comm);
        par::Mpi_Bcast(&EM3_BLK_MAX_Z,1,0,comm);

        EM3_OCTREE_MAX[0]=(double )(1u<<em3::EM3_MAXDEPTH);
        EM3_OCTREE_MAX[1]=(double )(1u<<em3::EM3_MAXDEPTH);
        EM3_OCTREE_MAX[2]=(double )(1u<<em3::EM3_MAXDEPTH);

        EM3_COMPD_MIN[0]=EM3_GRID_MIN_X;
        EM3_COMPD_MIN[1]=EM3_GRID_MIN_Y;
        EM3_COMPD_MIN[2]=EM3_GRID_MIN_Z;

        EM3_COMPD_MAX[0]=EM3_GRID_MAX_X;
        EM3_COMPD_MAX[1]=EM3_GRID_MAX_Y;
        EM3_COMPD_MAX[2]=EM3_GRID_MAX_Z;
        EM3_PADDING_WIDTH = EM3_ELE_ORDER>>1u;


        par::Mpi_Bcast(&KO_DISS_SIGMA, 1, 0, comm);
        par::Mpi_Bcast(&DISSIPATION_TYPE, 1, 0, comm);

        par::Mpi_Bcast(&EM3_NUM_REFINE_VARS,1,0,comm);
        par::Mpi_Bcast(&EM3_NUM_EVOL_VARS_VTU_OUTPUT,1,0,comm);
        par::Mpi_Bcast(&EM3_NUM_CONST_VARS_VTU_OUTPUT,1,0,comm);


        if(EM3_NUM_REFINE_VARS>EM3_NUM_VARS){std::cout<<"Error[parameter file]: Number of refine variables should be less than number of EM3_NUM_VARS"<<std::endl; exit(0);}
        if(EM3_NUM_EVOL_VARS_VTU_OUTPUT>EM3_NUM_VARS){std::cout<<"Error[parameter file]: Number of evolution VTU variables should be less than number of EM3_NUM_VARS"<<std::endl; exit(0);}

        par::Mpi_Bcast(EM3_REFINE_VARIABLE_INDICES,EM3_NUM_VARS,0,comm);
        par::Mpi_Bcast(EM3_VTU_OUTPUT_EVOL_INDICES,EM3_NUM_VARS,0,comm);
        par::Mpi_Bcast(EM3_VTU_OUTPUT_CONST_INDICES,EM3_CONSTRAINT_NUM_VARS,0,comm);

    }

    void dumpParamFile(std::ostream& sout, int root, MPI_Comm comm)
    {
        int rank,npes;
        MPI_Comm_rank(comm,&rank);
        MPI_Comm_size(comm,&npes);

        if(rank ==root)
        {
            sout<<"rank: "<<root<<" parameter values "<<std::endl;
            sout<<YLW<<"\tnpes :"<<npes<<NRM<<std::endl;
            sout<<YLW<<"\tEM3_ELE_ORDER :"<<em3::EM3_ELE_ORDER<<NRM<<std::endl;
            sout<<YLW<<"\tEM3_PADDING_WIDTH :"<<em3::EM3_PADDING_WIDTH<<NRM<<std::endl;
            sout<<YLW<<"\tEM3_DIM :"<<em3::EM3_DIM<<NRM<<std::endl;
            sout<<YLW<<"\tEM3_IO_OUTPUT_FREQ :"<<em3::EM3_IO_OUTPUT_FREQ<<NRM<<std::endl;
            sout<<YLW<<"\tEM3_REMESH_TEST_FREQ :"<<em3::EM3_REMESH_TEST_FREQ<<NRM<<std::endl;
            sout<<YLW<<"\tEM3_CHECKPT_FREQ :"<<em3::EM3_CHECKPT_FREQ<<NRM<<std::endl;
            sout<<YLW<<"\tEM3_RESTORE_SOLVER :"<<em3::EM3_RESTORE_SOLVER<<NRM<<std::endl;
            sout<<YLW<<"\tEM3_ENABLE_BLOCK_ADAPTIVITY :"<<em3::EM3_ENABLE_BLOCK_ADAPTIVITY<<NRM<<std::endl;
            sout<<YLW<<"\tEM3_VTU_FILE_PREFIX :"<<em3::EM3_VTU_FILE_PREFIX<<NRM<<std::endl;
            sout<<YLW<<"\tEM3_CHKPT_FILE_PREFIX :"<<em3::EM3_CHKPT_FILE_PREFIX<<NRM<<std::endl;
            sout<<YLW<<"\tEM3_PROFILE_FILE_PREFIX :"<<em3::EM3_PROFILE_FILE_PREFIX<<NRM<<std::endl;
            sout<<YLW<<"\tEM3_VTU_Z_SLICE_ONLY :"<<em3::EM3_VTU_Z_SLICE_ONLY<<NRM<<std::endl;
            sout<<YLW<<"\tEM3_IO_OUTPUT_GAP :"<<em3::EM3_IO_OUTPUT_GAP<<NRM<<std::endl;
            sout<<YLW<<"\tEM3_DENDRO_GRAIN_SZ :"<<em3::EM3_DENDRO_GRAIN_SZ<<NRM<<std::endl;
            sout<<YLW<<"\tEM3_ASYNC_COMM_K :"<<em3::EM3_ASYNC_COMM_K<<NRM<<std::endl;
            sout<<YLW<<"\tEM3_DENDRO_AMR_FAC :"<<em3::EM3_DENDRO_AMR_FAC<<NRM<<std::endl;
            sout<<YLW<<"\tEM3_CFL_FACTOR:"<<em3::EM3_CFL_FACTOR<<NRM<<std::endl;
            sout<<YLW<<"\tEM3_WAVELET_TOL :"<<em3::EM3_WAVELET_TOL<<NRM<<std::endl;
            sout<<YLW<<"\tEM3_LOAD_IMB_TOL :"<<em3::EM3_LOAD_IMB_TOL<<NRM<<std::endl;
            sout<<YLW<<"\tEM3_RK45_TIME_BEGIN :"<<em3::EM3_RK45_TIME_BEGIN<<NRM<<std::endl;
            sout<<YLW<<"\tEM3_RK45_TIME_END :"<<em3::EM3_RK45_TIME_END<<NRM<<std::endl;
            sout<<YLW<<"\tEM3_RK45_TIME_STEP_SIZE :"<<em3::EM3_RK45_TIME_STEP_SIZE<<NRM<<std::endl;
            sout<<YLW<<"\tEM3_RK45_DESIRED_TOL :"<<em3::EM3_RK45_DESIRED_TOL<<NRM<<std::endl;
            sout<<YLW<<"\tEM3_COMPD_MIN : ( :"<<em3::EM3_COMPD_MIN[0]<<" ,"<<em3::EM3_COMPD_MIN[1]<<","<<em3::EM3_COMPD_MIN[2]<<" )"<<NRM<<std::endl;
            sout<<YLW<<"\tEM3_COMPD_MAX : ( :"<<em3::EM3_COMPD_MAX[0]<<" ,"<<em3::EM3_COMPD_MAX[1]<<","<<em3::EM3_COMPD_MAX[2]<<" )"<<NRM<<std::endl;
            sout<<YLW<<"\tEM3_BLK_MIN : ( :"<<em3::EM3_BLK_MIN_X<<" ,"<<em3::EM3_BLK_MIN_Y<<","<<em3::EM3_BLK_MIN_Z<<" )"<<NRM<<std::endl;
            sout<<YLW<<"\tEM3_BLK_MAX : ( :"<<em3::EM3_BLK_MAX_X<<" ,"<<em3::EM3_BLK_MAX_Y<<","<<em3::EM3_BLK_MAX_Z<<" )"<<NRM<<std::endl;
            sout<<YLW<<"\tEM3_OCTREE_MIN : ( :"<<em3::EM3_OCTREE_MIN[0]<<" ,"<<em3::EM3_OCTREE_MIN[1]<<","<<em3::EM3_OCTREE_MIN[2]<<" )"<<NRM<<std::endl;
            sout<<YLW<<"\tEM3_OCTREE_MAX : ( :"<<em3::EM3_OCTREE_MAX[0]<<" ,"<<em3::EM3_OCTREE_MAX[1]<<","<<em3::EM3_OCTREE_MAX[2]<<" )"<<NRM<<std::endl;
            sout<<YLW<<"\tKO_DISS_SIGMA :"<<em3::KO_DISS_SIGMA<<NRM<<std::endl;
            sout<<YLW<<"\tEM3_ID_TYPE:"<<em3::EM3_ID_TYPE<<NRM<<std::endl;
            sout<<YLW<<"\tEM3_ID_AMP1:"<<em3::EM3_ID_AMP1<<NRM<<std::endl;
            sout<<YLW<<"\tEM3_ID_LAMBDA1:"<<em3::EM3_ID_LAMBDA1<<NRM<<std::endl;
            sout<<YLW<<"\tEM3_ID_AMP2:"<<em3::EM3_ID_AMP2<<NRM<<std::endl;
            //sout<<YLW<<"\tEM3_ID_DELTA1:"<<em3::EM3_ID_DELTA1<<NRM<<std::endl;
            //sout<<YLW<<"\tEM3_ID_DELTA2:"<<em3::EM3_ID_DELTA2<<NRM<<std::endl;
            //sout<<YLW<<"\tEM3_ID_XC1:"<<em3::EM3_ID_XC1<<NRM<<std::endl;
            //sout<<YLW<<"\tEM3_ID_YC1:"<<em3::EM3_ID_YC1<<NRM<<std::endl;
            //sout<<YLW<<"\tEM3_ID_ZC1:"<<em3::EM3_ID_ZC1<<NRM<<std::endl;
            //sout<<YLW<<"\tEM3_ID_XC2:"<<em3::EM3_ID_XC2<<NRM<<std::endl;
            //sout<<YLW<<"\tEM3_ID_YC2:"<<em3::EM3_ID_YC2<<NRM<<std::endl;
            //sout<<YLW<<"\tEM3_ID_ZC2:"<<em3::EM3_ID_ZC2<<NRM<<std::endl;
            //sout<<YLW<<"\tEM3_ID_EPSX1:"<<em3::EM3_ID_EPSX1<<NRM<<std::endl;
            //sout<<YLW<<"\tEM3_ID_EPSY1:"<<em3::EM3_ID_EPSY1<<NRM<<std::endl;
            //sout<<YLW<<"\tEM3_ID_EPSZ1:"<<em3::EM3_ID_EPSY1<<NRM<<std::endl;
            //sout<<YLW<<"\tEM3_ID_EPSX2:"<<em3::EM3_ID_EPSX2<<NRM<<std::endl;
            //sout<<YLW<<"\tEM3_ID_EPSY2:"<<em3::EM3_ID_EPSY2<<NRM<<std::endl;
            //sout<<YLW<<"\tEM3_ID_EPSZ2:"<<em3::EM3_ID_EPSY2<<NRM<<std::endl;
            //sout<<YLW<<"\tEM3_ID_R1:"<<em3::EM3_ID_R1<<NRM<<std::endl;
            //sout<<YLW<<"\tEM3_ID_R2:"<<em3::EM3_ID_R2<<NRM<<std::endl;
            //sout<<YLW<<"\tEM3_ID_NU1:"<<em3::EM3_ID_NU1<<NRM<<std::endl;
            //sout<<YLW<<"\tEM3_ID_NU2:"<<em3::EM3_ID_NU2<<NRM<<std::endl;
            //sout<<YLW<<"\tEM3_ID_OMEGA:"<<em3::EM3_ID_OMEGA<<NRM<<std::endl;

        
        

            //std::cout<<YLW<<"\tEM3_DIM :"<<em3::EM3_DIM<<NRM<<std::endl;
            sout<<YLW<<"\tEM3_MAXDEPTH :"<<em3::EM3_MAXDEPTH<<NRM<<std::endl;

            sout<<YLW<<"\tEM3_NUM_REFINE_VARS :"<<em3::EM3_NUM_REFINE_VARS<<NRM<<std::endl;
            sout<<YLW<<"\tEM3_REFINE_VARIABLE_INDICES :[";
            for(unsigned int i=0;i<em3::EM3_NUM_REFINE_VARS-1;i++)
                sout<<em3::EM3_REFINE_VARIABLE_INDICES[i]<<", ";
            sout<<em3::EM3_REFINE_VARIABLE_INDICES[em3::EM3_NUM_REFINE_VARS-1]<<"]"<<NRM<<std::endl;

            sout<<YLW<<"\tEM3_REFINEMENT_MODE :"<<em3::EM3_REFINEMENT_MODE<<NRM<<std::endl;

            sout<<YLW<<"\tEM3_NUM_EVOL_VARS_VTU_OUTPUT :"<<em3::EM3_NUM_EVOL_VARS_VTU_OUTPUT<<NRM<<std::endl;
            sout<<YLW<<"\tEM3_VTU_OUTPUT_EVOL_INDICES :[";
            for(unsigned int i=0;i<em3::EM3_NUM_EVOL_VARS_VTU_OUTPUT-1;i++)
                sout<<em3::EM3_VTU_OUTPUT_EVOL_INDICES[i]<<", ";
            sout<<em3::EM3_VTU_OUTPUT_EVOL_INDICES[em3::EM3_NUM_EVOL_VARS_VTU_OUTPUT-1]<<"]"<<NRM<<std::endl;


            #ifdef EM3_USE_4TH_ORDER_DERIVS
                sout<<"Using 4th order FD stencils. "<<std::endl;
            #endif

            #ifdef EM3_USE_6TH_ORDER_DERIVS
                sout<<"Using 6th order FD stencils. "<<std::endl;
            #endif

            #ifdef EM3_USE_8TH_ORDER_DERIVS
                sout<<"Using 8th order FD stencils. "<<std::endl;
            #endif

        }

        return ;
    }

    void initData(const double xx1,const double yy1,const double zz1, double *var)
    {
        const double x=GRIDX_TO_X(xx1);
        const double y=GRIDY_TO_Y(yy1);
        const double z=GRIDZ_TO_Z(zz1);

        const double amp1 = em3::EM3_ID_AMP1;
        const double lambda1 = em3::EM3_ID_LAMBDA1;

        double B0,B1,B2, E0,E1,E2; 
        double rho_e, J0,J1,J2;
        double r,Ephi_over_r,tmp_Ephiup ; 

        r = sqrt( x*x + y*y + z*z ) ; 
        tmp_Ephiup = - 8.0*amp1*lambda1*lambda1*exp(-lambda1*r*r) ; 
        E0 = - y * tmp_Ephiup ; 
        E1 =   x * tmp_Ephiup ; 
        E2 = 0.0 ; 
  
        B0 = 0.0 ;  
        B1 = 0.0 ;  
        B2 = 0.0 ;  
        
        J0 = 0.0 ;  
        J1 = 0.0 ;  
        J2 = 0.0 ;  
        
        rho_e = 0.0 ;  
        
        var[VAR::U_E0] = E0;
        var[VAR::U_E1] = E1;
        var[VAR::U_E2] = E2;
        var[VAR::U_B0] = B0;
        var[VAR::U_B1] = B1;
        var[VAR::U_B2] = B2;

    }


    void analyticalSol(const double xx1,const double yy1,const double zz1,const double t, double *var)
    {
        const double x=GRIDX_TO_X(xx1);
        const double y=GRIDY_TO_Y(yy1);
        const double z=GRIDZ_TO_Z(zz1);

        const double amp1 = em3::EM3_ID_AMP1;
        const double lambda1 = em3::EM3_ID_LAMBDA1;

        double B0,B1,B2, E0,E1,E2;
        double rho_e, J0,J1,J2;
        double r,tmp_Ephiup, tmp_Br_up, tmp_Btheta_up ;

        r = sqrt( x*x + y*y + z*z ) ;
        if ( r < 1.e-8 ) { r=1.e-8 ; }

        tmp_Br_up = - 2.0 * lambda1
                          * (   (t-r)*exp( - lambda1*(t-r)*(t-r) )
                              + (t+r)*exp( - lambda1*(t+r)*(t+r) ) ) / (r*r)
                    + (   exp( - lambda1*(t-r)*(t-r) )
                        - exp( - lambda1*(t+r)*(t+r) ) ) / (r*r*r)
                   ;
        tmp_Br_up *= 2.0 * amp1 ;  

        tmp_Btheta_up = - 2.0*lambda1
                             * (   exp( - lambda1*(t-r)*(t-r) )
                                 - exp( - lambda1*(t+r)*(t+r) ) ) / r
                        + 4.0*lambda1*lambda1
                             * (   (t-r)*(t-r)*exp( - lambda1*(t-r)*(t-r) )
                                 - (t+r)*(t+r)*exp( - lambda1*(t+r)*(t+r) ) ) 
                             / r
                        - 2.0*lambda1
                             * (   (t-r)*exp( - lambda1*(t-r)*(t-r) )
                                 + (t+r)*exp( - lambda1*(t+r)*(t+r) ) ) / (r*r)
                       + (   exp( - lambda1*(t-r)*(t-r) )
                           - exp( - lambda1*(t+r)*(t+r) ) ) / (r*r*r)
                   ;
        tmp_Btheta_up *= amp1 ;  

        tmp_Ephiup =   2.0*amp1*lambda1
                          * (   (t-r)*exp( - lambda1*(t-r)*(t-r) )
                              - (t+r)*exp( - lambda1*(t+r)*(t+r) ) ) / (r*r)
                     + 2.0*amp1*lambda1
                          * (   exp( - lambda1*(t-r)*(t-r) )
                              + exp( - lambda1*(t+r)*(t+r) ) ) / r
                     - 4.0*amp1*lambda1*lambda1
                          * (   (t-r)*(t-r)*exp( - lambda1*(t-r)*(t-r) )
                              + (t+r)*(t+r)*exp( - lambda1*(t+r)*(t+r) ) ) / r
                   ;

        E0 = - y * tmp_Ephiup / r ;
        E1 =   x * tmp_Ephiup / r ;
        E2 = 0.0 ;

        B0 = x*z * ( tmp_Br_up + tmp_Btheta_up ) / (r*r) ;
        B1 = y*z * ( tmp_Br_up + tmp_Btheta_up ) / (r*r) ;
        B2 = ( z*z * tmp_Br_up - (x*x+y*y) * tmp_Btheta_up ) / (r*r) ;

        J0 = 0.0 ;
        J1 = 0.0 ;
        J2 = 0.0 ;

        var[VAR::U_E0] = E0 ;
        var[VAR::U_E1] = E1 ;
        var[VAR::U_E2] = E2 ;
        var[VAR::U_B0] = B0 ;
        var[VAR::U_B1] = B1 ;
        var[VAR::U_B2] = B2 ;


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

    double computeWTol(double x,double y,double z,double* hx)
    {
       return em3::EM3_WAVELET_TOL;
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

    unsigned int getOctantWeight(const ot::TreeNode* pNode)
    {
        return (1u<<(3*pNode->getLevel()))*1;
    }

}// end of namespace em3





namespace em3
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
                outfile<<"partition tol : "<<em3::EM3_LOAD_IMB_TOL<<std::endl;
                outfile<<"wavelet tol : "<<em3::EM3_WAVELET_TOL<<std::endl;
                outfile<<"maxdepth : "<<em3::EM3_MAXDEPTH<<std::endl;

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


#ifdef EM3_PROFILE_HUMAN_READABLE
            if(!activeRank)
            {
                sprintf(fName,"%s_im.prof",filePrefix);
                outfile.open (fName,std::fstream::app);
                if(outfile.fail()) {std::cout<<fName<<" file open failed "<<std::endl; return ;}

                outfile<<"active npes : "<<activeNpes<<std::endl;
                outfile<<"global npes : "<<globalNpes<<std::endl;
                outfile<<"current step : "<<currentStep<<std::endl;
                outfile<<"partition tol : "<<em3::EM3_LOAD_IMB_TOL<<std::endl;
                outfile<<"wavelet tol : "<<em3::EM3_WAVELET_TOL<<std::endl;
                outfile<<"maxdepth : "<<em3::EM3_MAXDEPTH<<std::endl;

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
            if(!rank) outfile<<em3::EM3_LOAD_IMB_TOL<<"\t ";
            if(!rank) outfile<<em3::EM3_WAVELET_TOL<<"\t ";
            if(!rank) outfile<<em3::EM3_MAXDEPTH<<"\t ";

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
