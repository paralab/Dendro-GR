//
// Created by milinda on 7/26/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief Contains utility functions for MAXWELL simulation.
*/
//

#include "maxwellUtils.h"

namespace maxwell
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

            maxwell::MAXWELL_IO_OUTPUT_FREQ=parFile["MAXWELL_IO_OUTPUT_FREQ"];
            maxwell::MAXWELL_REMESH_TEST_FREQ=parFile["MAXWELL_REMESH_TEST_FREQ"];
            maxwell::MAXWELL_CHECKPT_FREQ=parFile["MAXWELL_CHECKPT_FREQ"];
            maxwell::MAXWELL_IO_OUTPUT_GAP=parFile["MAXWELL_IO_OUTPUT_GAP"];
            maxwell::MAXWELL_VTU_FILE_PREFIX=parFile["MAXWELL_VTU_FILE_PREFIX"].get<std::string>();
            maxwell::MAXWELL_CHKPT_FILE_PREFIX=parFile["MAXWELL_CHKPT_FILE_PREFIX"].get<std::string>();
            maxwell::MAXWELL_PROFILE_FILE_PREFIX=parFile["MAXWELL_PROFILE_FILE_PREFIX"].get<std::string>();
            maxwell::MAXWELL_RESTORE_SOLVER=parFile["MAXWELL_RESTORE_SOLVER"];
            maxwell::MAXWELL_ENABLE_BLOCK_ADAPTIVITY=parFile["MAXWELL_ENABLE_BLOCK_ADAPTIVITY"];
            maxwell::MAXWELL_BLK_MIN_X=parFile["MAXWELL_BLK_MIN_X"];
            maxwell::MAXWELL_BLK_MIN_Y=parFile["MAXWELL_BLK_MIN_Y"];
            maxwell::MAXWELL_BLK_MIN_Z=parFile["MAXWELL_BLK_MIN_Z"];
            maxwell::MAXWELL_BLK_MAX_X=parFile["MAXWELL_BLK_MAX_X"];
            maxwell::MAXWELL_BLK_MAX_Y=parFile["MAXWELL_BLK_MAX_Y"];
            maxwell::MAXWELL_BLK_MAX_Z=parFile["MAXWELL_BLK_MAX_Z"];
            maxwell::MAXWELL_DENDRO_GRAIN_SZ=parFile["MAXWELL_DENDRO_GRAIN_SZ"];
            maxwell::MAXWELL_ASYNC_COMM_K=parFile["MAXWELL_ASYNC_COMM_K"];
            maxwell::MAXWELL_DENDRO_AMR_FAC=parFile["MAXWELL_DENDRO_AMR_FAC"];
            maxwell::MAXWELL_LOAD_IMB_TOL=parFile["MAXWELL_LOAD_IMB_TOL"];
            maxwell::MAXWELL_RK45_TIME_BEGIN=parFile["MAXWELL_RK45_TIME_BEGIN"];
            maxwell::MAXWELL_RK45_TIME_END=parFile["MAXWELL_RK45_TIME_END"];
            maxwell::MAXWELL_RK45_TIME_STEP_SIZE=parFile["MAXWELL_RK45_TIME_STEP_SIZE"];
            maxwell::MAXWELL_RK45_DESIRED_TOL=parFile["MAXWELL_RK45_DESIRED_TOL"];
            maxwell::MAXWELL_DIM=parFile["MAXWELL_DIM"];
            maxwell::MAXWELL_MAXDEPTH=parFile["MAXWELL_MAXDEPTH"];
            maxwell::MAXWELL_GRID_MIN_X=parFile["MAXWELL_GRID_MIN_X"];
            maxwell::MAXWELL_GRID_MAX_X=parFile["MAXWELL_GRID_MAX_X"];
            maxwell::MAXWELL_GRID_MIN_Y=parFile["MAXWELL_GRID_MIN_Y"];
            maxwell::MAXWELL_GRID_MAX_Y=parFile["MAXWELL_GRID_MAX_Y"];
            maxwell::MAXWELL_GRID_MIN_Z=parFile["MAXWELL_GRID_MIN_Z"];
            maxwell::MAXWELL_GRID_MAX_Z=parFile["MAXWELL_GRID_MAX_Z"];
            maxwell::KO_DISS_SIGMA=parFile["KO_DISS_SIGMA"];

            maxwell::MAXWELL_ID_TYPE=parFile["MAXWELL_ID_TYPE"];

            maxwell::MAXWELL_WAVELET_TOL=parFile["MAXWELL_WAVELET_TOL"];

            maxwell::MAXWELL_NUM_REFINE_VARS=parFile["MAXWELL_NUM_REFINE_VARS"];
            for(unsigned int i=0;i<maxwell::MAXWELL_NUM_REFINE_VARS;i++)
                maxwell::MAXWELL_REFINE_VARIABLE_INDICES[i]=parFile["MAXWELL_REFINE_VARIABLE_INDICES"][i];

            maxwell::MAXWELL_NUM_EVOL_VARS_VTU_OUTPUT=parFile["MAXWELL_NUM_EVOL_VARS_VTU_OUTPUT"];

            for(unsigned int i=0;i<maxwell::MAXWELL_NUM_EVOL_VARS_VTU_OUTPUT;i++)
                maxwell::MAXWELL_VTU_OUTPUT_EVOL_INDICES[i]=parFile["MAXWELL_VTU_OUTPUT_EVOL_INDICES"][i];


            vtu_len=MAXWELL_VTU_FILE_PREFIX.size();
            chp_len=MAXWELL_CHKPT_FILE_PREFIX.size();
            prf_len=MAXWELL_PROFILE_FILE_PREFIX.size();

        }


        par::Mpi_Bcast(&MAXWELL_IO_OUTPUT_FREQ,1,0,comm);
        par::Mpi_Bcast(&MAXWELL_REMESH_TEST_FREQ,1,0,comm);
        par::Mpi_Bcast(&MAXWELL_CHECKPT_FREQ,1,0,comm);
        par::Mpi_Bcast(&MAXWELL_IO_OUTPUT_GAP,1,0,comm);

        par::Mpi_Bcast(&vtu_len,1,0,comm);
        par::Mpi_Bcast(&chp_len,1,0,comm);
        par::Mpi_Bcast(&prf_len,1,0,comm);

        par::Mpi_Bcast(&MAXWELL_DENDRO_GRAIN_SZ,1,0,comm);
        par::Mpi_Bcast(&MAXWELL_DENDRO_AMR_FAC,1,0,comm);
        par::Mpi_Bcast(&MAXWELL_ASYNC_COMM_K,1,0,comm);

        char vtu_name[vtu_len+1];
        char chp_name[chp_len+1];
        char prf_name[prf_len+1];


        if(!rank)
        {
           for(unsigned int k=0;k<vtu_len;k++)
               vtu_name[k]=MAXWELL_VTU_FILE_PREFIX[k];

            for(unsigned int k=0;k<chp_len;k++)
                chp_name[k]=MAXWELL_CHKPT_FILE_PREFIX[k];

            for(unsigned int k=0;k<prf_len;k++)
                prf_name[k]=MAXWELL_PROFILE_FILE_PREFIX[k];

            vtu_name[vtu_len]='\0';
            chp_name[chp_len]='\0';
            prf_name[prf_len]='\0';

        }


        MPI_Bcast(vtu_name,vtu_len+1,MPI_CHAR,0,comm);
        MPI_Bcast(chp_name,chp_len+1,MPI_CHAR,0,comm);
        MPI_Bcast(prf_name,prf_len+1,MPI_CHAR,0,comm);

        MAXWELL_VTU_FILE_PREFIX=std::string(vtu_name);
        MAXWELL_CHKPT_FILE_PREFIX=std::string(chp_name);
        MAXWELL_PROFILE_FILE_PREFIX=std::string(prf_name);


        par::Mpi_Bcast(&MAXWELL_RESTORE_SOLVER,1,0,comm);
        par::Mpi_Bcast(&MAXWELL_ENABLE_BLOCK_ADAPTIVITY,1,0,comm);

        par::Mpi_Bcast(&MAXWELL_WAVELET_TOL,1,0,comm);

        par::Mpi_Bcast(&MAXWELL_LOAD_IMB_TOL,1,0,comm);
        par::Mpi_Bcast(&MAXWELL_RK45_TIME_BEGIN,1,0,comm);
        par::Mpi_Bcast(&MAXWELL_RK45_TIME_END,1,0,comm);
        par::Mpi_Bcast(&MAXWELL_RK45_TIME_STEP_SIZE,1,0,comm);
        par::Mpi_Bcast(&MAXWELL_RK45_DESIRED_TOL,1,0,comm);
        par::Mpi_Bcast(&MAXWELL_DIM,1,0,comm);
        par::Mpi_Bcast(&MAXWELL_MAXDEPTH,1,0,comm);

        par::Mpi_Bcast(&MAXWELL_ID_TYPE,1,0,comm);

        par::Mpi_Bcast(&MAXWELL_GRID_MIN_X,1,0,comm);
        par::Mpi_Bcast(&MAXWELL_GRID_MAX_X,1,0,comm);
        par::Mpi_Bcast(&MAXWELL_GRID_MIN_Y,1,0,comm);
        par::Mpi_Bcast(&MAXWELL_GRID_MAX_Y,1,0,comm);
        par::Mpi_Bcast(&MAXWELL_GRID_MIN_Z,1,0,comm);
        par::Mpi_Bcast(&MAXWELL_GRID_MAX_Z,1,0,comm);

        par::Mpi_Bcast(&MAXWELL_BLK_MIN_X,1,0,comm);
        par::Mpi_Bcast(&MAXWELL_BLK_MIN_Y,1,0,comm);
        par::Mpi_Bcast(&MAXWELL_BLK_MIN_Z,1,0,comm);

        par::Mpi_Bcast(&MAXWELL_BLK_MAX_X,1,0,comm);
        par::Mpi_Bcast(&MAXWELL_BLK_MAX_Y,1,0,comm);
        par::Mpi_Bcast(&MAXWELL_BLK_MAX_Z,1,0,comm);

        MAXWELL_OCTREE_MAX[0]=(double )(1u<<maxwell::MAXWELL_MAXDEPTH);
        MAXWELL_OCTREE_MAX[1]=(double )(1u<<maxwell::MAXWELL_MAXDEPTH);
        MAXWELL_OCTREE_MAX[2]=(double )(1u<<maxwell::MAXWELL_MAXDEPTH);

        MAXWELL_COMPD_MIN[0]=MAXWELL_GRID_MIN_X;
        MAXWELL_COMPD_MIN[1]=MAXWELL_GRID_MIN_Y;
        MAXWELL_COMPD_MIN[2]=MAXWELL_GRID_MIN_Z;

        MAXWELL_COMPD_MAX[0]=MAXWELL_GRID_MAX_X;
        MAXWELL_COMPD_MAX[1]=MAXWELL_GRID_MAX_Y;
        MAXWELL_COMPD_MAX[2]=MAXWELL_GRID_MAX_Z;


        par::Mpi_Bcast(&KO_DISS_SIGMA, 1, 0, comm);

        par::Mpi_Bcast(&MAXWELL_NUM_REFINE_VARS,1,0,comm);
        par::Mpi_Bcast(&MAXWELL_NUM_EVOL_VARS_VTU_OUTPUT,1,0,comm);


        if(MAXWELL_NUM_REFINE_VARS>MAXWELL_NUM_VARS){std::cout<<"Error[parameter file]: Number of refine variables should be less than number of MAXWELL_NUM_VARS"<<std::endl; exit(0);}
        if(MAXWELL_NUM_EVOL_VARS_VTU_OUTPUT>MAXWELL_NUM_VARS){std::cout<<"Error[parameter file]: Number of evolution VTU variables should be less than number of MAXWELL_NUM_VARS"<<std::endl; exit(0);}

        par::Mpi_Bcast(MAXWELL_REFINE_VARIABLE_INDICES,MAXWELL_NUM_VARS,0,comm);
        par::Mpi_Bcast(MAXWELL_VTU_OUTPUT_EVOL_INDICES,MAXWELL_NUM_VARS,0,comm);

    }



    void initData(const double xx1,const double yy1,const double zz1, double *var)
    {

        const double x=GRIDX_TO_X(xx1);
        const double y=GRIDY_TO_Y(yy1);
        const double z=GRIDZ_TO_Z(zz1);
        const double r=sqrt(x*x + y*y + z*z);

        //std::cout<<"initData: "<<x<<", "<<y<<", "<<z<<std::endl;

        if (maxwell::MAXWELL_ID_TYPE == 0) {
          var[VAR::U_AX] = 0;
          var[VAR::U_AY] = 0;
          var[VAR::U_AZ] = 0;
	      var[VAR::U_EX] = -8.0*y*exp(-r*r);
	      var[VAR::U_EY] = 8.0*x*exp(-r*r);
	      var[VAR::U_EZ] = 0;
	      var[VAR::U_GAM] = 0;
	      var[VAR::U_PSI] = 0;
        }

	else if (maxwell::MAXWELL_ID_TYPE == 1) {
          var[VAR::U_AX] = 0;
          var[VAR::U_AY] = 0;
          var[VAR::U_AZ] = 0;
	      var[VAR::U_EX] = .001*sin(3.0*x)+.001*sin(5.0*x)+.001*sin(11.0*x);
	      var[VAR::U_EY] = 0;
	      var[VAR::U_EZ] = 0;
	      var[VAR::U_GAM] = 0;
	      var[VAR::U_PSI] = 0;
        }

    }

    void solToroidalDipole(const double t, const double xx1, const double yy1, const double zz1, double *var)
    {
        const double x=GRIDX_TO_X(xx1);
        const double y=GRIDY_TO_Y(yy1);
        const double z=GRIDZ_TO_Z(zz1);
        const double r=sqrt(x*x + y*y + z*z);

        double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t16;
        double t20,t21,t22,t24,t30,t31,t39,t44,t47,t48,t51,t54;
        double t55,t63,t66,t69,t78,t81;
        // Ax
        t1 = x * x;
        t2 = y * y;
        t3 = z * z;
        t4 = t1 + t2 + t3 + 1.0e-11;
        t5 = sqrt(t4);
        t6 = 0.1e1 / t5;
        t8 = t - t5;
        t9 = t8 * t8;
        t10 = exp(-t9);
        t11 = t + t5;
        t12 = t11 * t11;
        t13 = exp(-t12);
        var[VAR::U_AX] = -(0.1e1 / t4 * (t10 - t13) - 0.2e1 * t6 * (t10 * t8 + t13 * t11)) * t6 * y;

        // Ay
        t6 = 0.1e1 / t5;
        t8 = t - t5;
        t9 = t8 * t8;
        t10 = exp(-t9);
        t11 = t + t5;
        t12 = t11 * t11;
        t13 = exp(-t12);
        var[VAR::U_AY] = (0.1e1 / t4 * (t10 - t13) - 0.2e1 * t6 * (t10 * t8 + t13 * t11)) * t6 * x;

        // Ex
        t7 = (t + t5)*(t + t5);
        t8 = exp(-t7);
        t9 = t * t8;
        t11 = t * t;
        t12 = t11 * t8;
        t20 = (-t + t5) * (-t + t5);
        t21 = exp(-t20);
        t22 = t * t21;
        t24 = t11 * t21;
        t31 = t1 * t8;
        t39 = t1 * t21;
        t44 = t5 * t9 + 0.2e1 * t1 * t12 + 0.2e1 * t2 * t12 + 0.2e1 * t3 * t12 - t5 * t22 + 0.2e1 * t1 * t24 + 0.2e1 * t2 * t24 + 0.2e1 * t3 * t24 + 0.4e1 * t2 * t31 + 0.4e1 * t3 * t31 + 0.4e1 * t3 * t2 * t8 + 0.4e1 * t2 * t39 + 0.4e1 * t3 * t39;
        t48 = t1 * t1;
        t51 = t2 * t2;
        t54 = t3 * t3;
        t63 = t1 * t5;
        t66 = t2 * t5;
        t69 = t3 * t5;
        t78 = 0.4e1 * t3 * t2 * t21 + 0.2e1 * t48 * t8 + 0.2e1 * t51 * t8 + 0.2e1 * t54 * t8 + 0.2e1 * t48 * t21 + 0.2e1 * t51 * t21 + 0.2e1 * t54 * t21 + 0.4e1 * t63 * t9 + 0.4e1 * t66 * t9 + 0.4e1 * t69 * t9 - 0.4e1 * t63 * t22 - 0.4e1 * t66 * t22 - 0.4e1 * t69 * t22;
        t81 = t4 * t4;
        var[VAR::U_EX] = -0.2e1 / t81 * (t44 + t78) * y;

        // Ey
        t7 = (t + t5) * (t + t5);
        t8 = exp(-t7);
        t9 = t * t8;
        t10 = t1 * t5;
        t13 = t2 * t5;
        t16 = t3 * t5;
        t20 = (-t + t5) * (-t + t5);
        t21 = exp(-t20);
        t22 = t * t21;
        t30 = t * t;
        t31 = t30 * t8;
        t39 = t30 * t21;
        t44 = 0.4e1 * t10 * t9 + 0.4e1 * t13 * t9 + 0.4e1 * t16 * t9 - 0.4e1 * t10 * t22 - 0.4e1 * t13 * t22 - 0.4e1 * t16 * t22 + t5 * t9 + 0.2e1 * t1 * t31 + 0.2e1 * t2 * t31 + 0.2e1 * t3 * t31 - t5 * t22 + 0.2e1 * t1 * t39 + 0.2e1 * t2 * t39;
        t47 = t1 * t8;
        t55 = t1 * t21;
        t63 = t1 * t1;
        t66 = t2 * t2;
        t69 = t3 * t3;
        t78 = 0.2e1 * t3 * t39 + 0.4e1 * t2 * t47 + 0.4e1 * t3 * t47 + 0.4e1 * t3 * t2 * t8 + 0.4e1 * t2 * t55 + 0.4e1 * t3 * t55 + 0.4e1 * t3 * t2 * t21 + 0.2e1 * t63 * t8 + 0.2e1 * t66 * t8 + 0.2e1 * t69 * t8 + 0.2e1 * t63 * t21 + 0.2e1 * t66 * t21 + 0.2e1 * t69 * t21;
        t81 = t4 * t4;
        var[VAR::U_EY] = 0.2e1 / t81 * (t44 + t78) * x;

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
       return maxwell::MAXWELL_WAVELET_TOL;
    }


}// end of namespace maxwell





namespace maxwell
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
                outfile<<"partition tol : "<<maxwell::MAXWELL_LOAD_IMB_TOL<<std::endl;
                outfile<<"wavelet tol : "<<maxwell::MAXWELL_WAVELET_TOL<<std::endl;
                outfile<<"maxdepth : "<<maxwell::MAXWELL_MAXDEPTH<<std::endl;

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


#ifdef MAXWELL_PROFILE_HUMAN_READABLE
            if(!activeRank)
            {
                sprintf(fName,"%s_im.prof",filePrefix);
                outfile.open (fName,std::fstream::app);
                if(outfile.fail()) {std::cout<<fName<<" file open failed "<<std::endl; return ;}

                outfile<<"active npes : "<<activeNpes<<std::endl;
                outfile<<"global npes : "<<globalNpes<<std::endl;
                outfile<<"current step : "<<currentStep<<std::endl;
                outfile<<"partition tol : "<<maxwell::MAXWELL_LOAD_IMB_TOL<<std::endl;
                outfile<<"wavelet tol : "<<maxwell::MAXWELL_WAVELET_TOL<<std::endl;
                outfile<<"maxdepth : "<<maxwell::MAXWELL_MAXDEPTH<<std::endl;

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
            if(!rank) outfile<<maxwell::MAXWELL_LOAD_IMB_TOL<<"\t ";
            if(!rank) outfile<<maxwell::MAXWELL_WAVELET_TOL<<"\t ";
            if(!rank) outfile<<maxwell::MAXWELL_MAXDEPTH<<"\t ";

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
