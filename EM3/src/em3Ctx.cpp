/**
 * @file em3Ctx.cpp
 * @author Milinda Fernando (milinda@cs.utah.edu)
 * @brief EM3 Ctx class. 
 * @version 0.1
 * @date 2020-07-18
 * 
 * @copyright Copyright (c) 2020
 * 
 */


#include "em3Ctx.h"
namespace em3
{
    EM3Ctx::EM3Ctx(ot::Mesh* pMesh) : Ctx()
    {
        m_uiMesh = pMesh;
        // variable allocation for evolution variables
        m_uiEVar = this->create_vec(ts::CTXVType::EVOLUTION,true,false,false,EM3_NUM_VARS);
        m_uiCVar = this->create_vec(ts::CTXVType::CONSTRAINT,true,false,false,EM3_CONSTRAINT_NUM_VARS);
        m_uiPVar = this->create_vec(ts::CTXVType::PRIMITIVE,true,false,false,EM3_NUM_VARS); // used as to compare the analytical solution with the evolved one. 

        m_uiEUnzip[0] = this->create_vec(ts::CTXVType::EVOLUTION,false,true,false,EM3_NUM_VARS);
        m_uiEUnzip[1] = this->create_vec(ts::CTXVType::EVOLUTION,false,true,false,EM3_NUM_VARS);

        m_uiCUnzip[0] = this->create_vec(ts::CTXVType::CONSTRAINT,false,true,false,EM3_CONSTRAINT_NUM_VARS);
        m_uiCUnzip[1] = this->create_vec(ts::CTXVType::CONSTRAINT,false,true,false,EM3_CONSTRAINT_NUM_VARS);

        m_uiTinfo._m_uiStep=0;
        m_uiTinfo._m_uiT = 0;
        m_uiTinfo._m_uiTb = EM3_RK45_TIME_BEGIN;
        m_uiTinfo._m_uiTe = EM3_RK45_TIME_END;
        m_uiTinfo._m_uiTh = EM3_RK45_TIME_STEP_SIZE;     

        m_uiElementOrder = EM3_ELE_ORDER;

        m_uiMinPt = Point(EM3_GRID_MIN_X,EM3_GRID_MIN_Y,EM3_GRID_MIN_Z);
        m_uiMaxPt = Point(EM3_GRID_MAX_X,EM3_GRID_MAX_Y,EM3_GRID_MAX_Z);

        return;

    }

    EM3Ctx::~EM3Ctx()
    {
        this->destroy_vec(m_uiEVar);
        this->destroy_vec(m_uiCVar);
        this->destroy_vec(m_uiPVar);

        this->destroy_vec(m_uiEUnzip[0]);
        this->destroy_vec(m_uiEUnzip[1]);

        this->destroy_vec(m_uiCUnzip[0]);
        this->destroy_vec(m_uiCUnzip[1]);


    }

    int EM3Ctx::initialize()
    {

        if(EM3_RESTORE_SOLVER)
        {
            this->restore_checkpt();
            return 0; 
        }

        unsigned int nodeLookUp_CG;
        unsigned int nodeLookUp_DG;
        DendroScalar x,y,z,len;
        const ot::TreeNode * pNodes=&(*(m_uiMesh->getAllElements().begin()));
        unsigned int ownerID,ii_x,jj_y,kk_z;
        unsigned int eleOrder=m_uiMesh->getElementOrder();
        const unsigned int * e2n_cg=&(*(m_uiMesh->getE2NMapping().begin()));
        const unsigned int * e2n_dg=&(*(m_uiMesh->getE2NMapping_DG().begin()));
        const unsigned int nPe=m_uiMesh->getNumNodesPerElement();
        const unsigned int nodeLocalBegin=m_uiMesh->getNodeLocalBegin();
        const unsigned int nodeLocalEnd=m_uiMesh->getNodeLocalEnd();

        const unsigned int numEVars = EM3_NUM_VARS;
        DendroScalar* var    = new DendroScalar [numEVars];
        DendroScalar** zipIn = new DendroScalar*[numEVars];
        m_uiEVar.Get2DArray(zipIn,true);

        DendroScalar mp, mm, mp_adm, mm_adm, E, J1, J2, J3;
        
        for(unsigned int elem=m_uiMesh->getElementLocalBegin(); elem<m_uiMesh->getElementLocalEnd(); elem++)
        {
            for(unsigned int k=0; k<(eleOrder+1); k++)
                for(unsigned int j=0; j<(eleOrder+1); j++ )
                    for(unsigned int i=0; i<(eleOrder+1); i++)
                    {
                        nodeLookUp_CG=e2n_cg[elem*nPe+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i];
                        if(nodeLookUp_CG>=nodeLocalBegin && nodeLookUp_CG<nodeLocalEnd)
                        {
                            nodeLookUp_DG=e2n_dg[elem*nPe+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i];
                            m_uiMesh->dg2eijk(nodeLookUp_DG,ownerID,ii_x,jj_y,kk_z);
                            len=(double)(1u<<(m_uiMaxDepth-pNodes[ownerID].getLevel()));
                            x=pNodes[ownerID].getX()+ ii_x*(len/(eleOrder));
                            y=pNodes[ownerID].getY()+ jj_y*(len/(eleOrder));
                            z=pNodes[ownerID].getZ()+ kk_z*(len/(eleOrder));

                            initData(x, y, z , var);

                            for(unsigned int v=0; v<numEVars; v++)
                                zipIn[v][nodeLookUp_CG]=var[v];
                            
                        }

                    }
        }
        
        delete [] var;
        delete [] zipIn;







        return 0;

    }

    int EM3Ctx::finalize()
    {
        return 0;    
    }

    int EM3Ctx::rhs(DVec* in , DVec* out, unsigned int sz , DendroScalar time)
    {

        // all the variables should be packed together. 
        assert(sz==1);

        DendroScalar ** sVar;
        in[0].Get2DArray(sVar,false);

        this->unzip(in[0],m_uiEUnzip[0] , EM3_ASYNC_COMM_K);
        
        DendroScalar ** unzipIn;
        DendroScalar **  unzipOut; 
        
        m_uiEUnzip[0].Get2DArray(unzipIn , false);
        m_uiEUnzip[1].Get2DArray(unzipOut ,false);

        const ot::Block* blkList = m_uiMesh->getLocalBlockList().data();
        const unsigned int numBlocks = m_uiMesh->getLocalBlockList().size();
        
        unsigned int offset;
        double ptmin[3], ptmax[3];
        unsigned int lsz[3];
        unsigned int bflag;
        double dx,dy,dz;
        
        const Point pt_min(em3::EM3_COMPD_MIN[0],em3::EM3_COMPD_MIN[1],em3::EM3_COMPD_MIN[2]);
        const Point pt_max(em3::EM3_COMPD_MAX[0],em3::EM3_COMPD_MAX[1],em3::EM3_COMPD_MAX[2]);
        const unsigned int PW=em3::EM3_PADDING_WIDTH;

        for(unsigned int blk =0; blk < numBlocks; blk++)
        {
            offset=blkList[blk].getOffset();
            lsz[0]=blkList[blk].getAllocationSzX();
            lsz[1]=blkList[blk].getAllocationSzY();
            lsz[2]=blkList[blk].getAllocationSzZ();

            bflag=blkList[blk].getBlkNodeFlag();

            dx=blkList[blk].computeDx(pt_min,pt_max);
            dy=blkList[blk].computeDy(pt_min,pt_max);
            dz=blkList[blk].computeDz(pt_min,pt_max);

            ptmin[0]=GRIDX_TO_X(blkList[blk].getBlockNode().minX())-PW*dx;
            ptmin[1]=GRIDY_TO_Y(blkList[blk].getBlockNode().minY())-PW*dy;
            ptmin[2]=GRIDZ_TO_Z(blkList[blk].getBlockNode().minZ())-PW*dz;

            ptmax[0]=GRIDX_TO_X(blkList[blk].getBlockNode().maxX())+PW*dx;
            ptmax[1]=GRIDY_TO_Y(blkList[blk].getBlockNode().maxY())+PW*dy;
            ptmax[2]=GRIDZ_TO_Z(blkList[blk].getBlockNode().maxZ())+PW*dz;

            #ifdef __PROFILE_CTX__
                m_uiCtxpt[ts::CTXPROFILE::RHS].start();
            #endif

            em3rhs(unzipOut,(const DendroScalar**)unzipIn, offset, ptmin, ptmax, lsz, bflag);
            
            // std::cout<<":::\n";
            // for(unsigned int n=0; n < lsz[0]*lsz[1]*lsz[2]; n++)
            //     std::cout<<"unzipOut: "<<n<<" val: "<<unzipOut[0][n]<<std::endl;

            #ifdef __PROFILE_CTX__
                m_uiCtxpt[ts::CTXPROFILE::RHS].stop();
            #endif


        }

        this->zip(m_uiEUnzip[1], out[0], EM3_ASYNC_COMM_K);

        // for(unsigned int v =0; v < in[0].GetDof(); v++)
        // {
        //     for(unsigned int n=0; n < m_uiMesh->getDegOfFreedom(); n++)
        //     {
        //         std::cout<<" cg::v: "<<v<<" n: "<<n<< " val: "<<(out[0].GetVecArray()+v*m_uiMesh->getDegOfFreedom())[n]<<std::endl;;
        //     }
        // }

        delete [] unzipIn;
        delete [] unzipOut;
        delete [] sVar;

        

        return 0;

    }

    int EM3Ctx::rhs_blk(const DendroScalar* in, DendroScalar* out, unsigned int dof, unsigned int local_blk_id, DendroScalar  blk_time)  
    {

        #ifdef __PROFILE_CTX__
            m_uiCtxpt[ts::CTXPROFILE::RHS_BLK].start();
        #endif

        // all the variables should be packed together. 
        const ot::Block* blkList = m_uiMesh->getLocalBlockList().data();
        const unsigned int numBlocks = m_uiMesh->getLocalBlockList().size();
        assert(local_blk_id < numBlocks);

        const unsigned int blk = local_blk_id;
        DendroScalar **  unzipIn   = new DendroScalar*[dof];
        DendroScalar **  unzipOut  = new DendroScalar*[dof]; 
        
        unsigned int offset;
        double ptmin[3], ptmax[3];
        unsigned int lsz[3];
        unsigned int bflag;
        double dx,dy,dz;
        
        const Point pt_min(em3::EM3_COMPD_MIN[0],em3::EM3_COMPD_MIN[1],em3::EM3_COMPD_MIN[2]);
        const Point pt_max(em3::EM3_COMPD_MAX[0],em3::EM3_COMPD_MAX[1],em3::EM3_COMPD_MAX[2]);
        const unsigned int PW=em3::EM3_PADDING_WIDTH;


        offset=blkList[blk].getOffset();
        lsz[0]=blkList[blk].getAllocationSzX();
        lsz[1]=blkList[blk].getAllocationSzY();
        lsz[2]=blkList[blk].getAllocationSzZ();

        const unsigned int NN  = lsz[0] * lsz[1] * lsz[2];

        for(unsigned int v =0; v < dof; v++)
        {
            unzipIn[v] = (DendroScalar*) (in + v*NN);
            unzipOut[v] = (DendroScalar*) (out + v*NN);
        }


        bflag=blkList[blk].getBlkNodeFlag();

        dx=blkList[blk].computeDx(pt_min,pt_max);
        dy=blkList[blk].computeDy(pt_min,pt_max);
        dz=blkList[blk].computeDz(pt_min,pt_max);

        ptmin[0]=GRIDX_TO_X(blkList[blk].getBlockNode().minX())-PW*dx;
        ptmin[1]=GRIDY_TO_Y(blkList[blk].getBlockNode().minY())-PW*dy;
        ptmin[2]=GRIDZ_TO_Z(blkList[blk].getBlockNode().minZ())-PW*dz;

        ptmax[0]=GRIDX_TO_X(blkList[blk].getBlockNode().maxX())+PW*dx;
        ptmax[1]=GRIDY_TO_Y(blkList[blk].getBlockNode().maxY())+PW*dy;
        ptmax[2]=GRIDZ_TO_Z(blkList[blk].getBlockNode().maxZ())+PW*dz;

        // note that the offset zero is important since it is the block vector. 
        em3rhs(unzipOut,(const DendroScalar**)unzipIn, 0, ptmin, ptmax, lsz, bflag);

        // for(unsigned int v =0; v < dof; v++)
        // {
        //     unsigned int nid=0; 
        //     for(unsigned int k=3; k < lsz[2]-3; k++)
        //     for(unsigned int j=3; j < lsz[1]-3; j++)
        //     for(unsigned int i=3; i < lsz[0]-3; i++,nid++)
        //     {
        //         std::cout<<" blk::v: "<<v<<" n: "<<nid<< " val: "<<unzipOut[v][k*lsz[1]*lsz[0] +j*lsz[0] +i ]<<std::endl;;
        //     }
        // }


        delete [] unzipIn;
        delete [] unzipOut;

        #ifdef __PROFILE_CTX__
            m_uiCtxpt[ts::CTXPROFILE::RHS_BLK].stop();
        #endif
        
        return 0;

    }

    int EM3Ctx::write_vtu()
    {

        if(!m_uiMesh->isActive())
            return 0;
        
        DendroScalar * consVar[EM3_CONSTRAINT_NUM_VARS];
        DendroScalar * evolVar[EM3_NUM_VARS];
        DendroScalar * primVar[EM3_NUM_VARS];

        m_uiEVar.Get2DVec(evolVar);
        m_uiPVar.Get2DVec(primVar);
        m_uiCVar.Get2DVec(consVar);

        
        this->compute_primitives();   // pvars are exhanged when computing the prims. 
        this->compute_constraints(); // evars and cvars are commuicated when computing the constraints. 
        

        std::vector<std::string> pDataNames;
        const unsigned int numConstVars = EM3_NUM_CONST_VARS_VTU_OUTPUT;
        const unsigned int numEvolVars  = EM3_NUM_EVOL_VARS_VTU_OUTPUT;
        const unsigned int numPrimVars  = EM3_NUM_EVOL_VARS_VTU_OUTPUT;

        
        double *pData[(numConstVars+numEvolVars + numPrimVars)];

        for(unsigned int i=0; i<numConstVars; i++)
        {
            pDataNames.push_back(std::string(EM3_CONSTRAINT_VAR_NAMES[EM3_VTU_OUTPUT_CONST_INDICES[i]]));
            pData[i]=consVar[EM3_VTU_OUTPUT_CONST_INDICES[i]];
        }

        for(unsigned int i=0; i<numEvolVars; i++)
        {
            pDataNames.push_back(std::string(EM3_VAR_NAMES[EM3_VTU_OUTPUT_EVOL_INDICES[i]]));
            pData[numConstVars + i] = evolVar[EM3_VTU_OUTPUT_EVOL_INDICES[i]];
        }

        for(unsigned int i=0; i<numPrimVars; i++)
        {
            pDataNames.push_back(std::string("diff_") + std::string(EM3_VAR_NAMES[EM3_VTU_OUTPUT_EVOL_INDICES[i]]));
            pData[numConstVars + numEvolVars + i] = primVar[EM3_VTU_OUTPUT_EVOL_INDICES[i]];
        }

        std::vector<char*> pDataNames_char;
        pDataNames_char.reserve(pDataNames.size());

        for(unsigned int  i = 0; i < pDataNames.size(); i++)
            pDataNames_char.push_back(const_cast<char*>(pDataNames[i].c_str()));

        const char * fDataNames[]= {"Time","Cycle"};
        const double fData[]= {m_uiTinfo._m_uiT,(double)m_uiTinfo._m_uiStep};

        char fPrefix[256];
        sprintf(fPrefix,"%s_%d",EM3_VTU_FILE_PREFIX.c_str(),m_uiTinfo._m_uiStep);

        io::vtk::mesh2vtuFine(m_uiMesh,fPrefix,2,fDataNames,fData,(numConstVars + numEvolVars + numPrimVars),(const char **)&pDataNames_char[0],(const double **)pData);
        return 0;

    }

    int EM3Ctx::write_checkpt()
    {   
        if(m_uiMesh->isActive())
        {
            unsigned int cpIndex;
            (m_uiTinfo._m_uiStep %(2*EM3_CHECKPT_FREQ)==0) ? cpIndex=0 : cpIndex=1; // to support alternate file writing.
            unsigned int rank=m_uiMesh->getMPIRank();
            unsigned int npes=m_uiMesh->getMPICommSize();

            char fName[256];
            const ot::TreeNode * pNodes=&(*(m_uiMesh->getAllElements().begin()+m_uiMesh->getElementLocalBegin()));
            sprintf(fName,"%s_octree_%d_%d.oct",EM3_CHKPT_FILE_PREFIX.c_str(),cpIndex,rank);
            io::checkpoint::writeOctToFile(fName,pNodes,m_uiMesh->getNumLocalMeshElements());

            unsigned int numVars   = EM3_NUM_VARS;
            const char ** varNames = EM3_VAR_NAMES;

            const unsigned int dof = m_uiEVar.GetDof();
            DendroScalar* eVar[dof];
            m_uiEVar.Get2DVec(eVar);

            sprintf(fName,"%s_%d_%d.var",EM3_CHKPT_FILE_PREFIX.c_str(),cpIndex,rank);
            io::checkpoint::writeVecToFile(fName,m_uiMesh,(const double **)eVar,em3::EM3_NUM_VARS);


            if(!rank)
            {
                sprintf(fName,"%s_step_%d.cp",EM3_CHKPT_FILE_PREFIX.c_str(),cpIndex);
                //std::cout<<"writing : "<<fName<<std::endl;
                std::ofstream outfile(fName);
                if(!outfile) {std::cout<<fName<<" file open failed "<<std::endl; return 0;}

                json checkPoint;

                checkPoint["DENDRO_TS_TIME_BEGIN"]    = m_uiTinfo._m_uiTb;
                checkPoint["DENDRO_TS_TIME_END"]      = m_uiTinfo._m_uiTe;
                checkPoint["DENDRO_TS_ELEMENT_ORDER"] = m_uiElementOrder;

                checkPoint["DENDRO_TS_TIME_CURRENT"]    = m_uiTinfo._m_uiT;
                checkPoint["DENDRO_TS_STEP_CURRENT"]    = m_uiTinfo._m_uiStep;
                checkPoint["DENDRO_TS_TIME_STEP_SIZE"]  = m_uiTinfo._m_uiTh;
                checkPoint["DENDRO_TS_LAST_IO_TIME"]    = m_uiTinfo._m_uiT;

                checkPoint["DENDRO_TS_WAVELET_TOLERANCE"]  = EM3_WAVELET_TOL;
                checkPoint["DENDRO_TS_LOAD_IMB_TOLERANCE"] = EM3_LOAD_IMB_TOL;
                checkPoint["DENDRO_TS_NUM_VARS"]=numVars; // number of variables to restore.
                checkPoint["DENDRO_TS_ACTIVE_COMM_SZ"]=m_uiMesh->getMPICommSize(); // (note that rank 0 is always active).
                
                outfile<<std::setw(4)<<checkPoint<<std::endl;
                outfile.close();

            }

        }

        return 0;
    }

    int EM3Ctx::restore_checkpt()
    {
        unsigned int numVars=0;
        std::vector<ot::TreeNode> octree;
        json checkPoint;

        int rank;
        int npes;
        MPI_Comm comm=m_uiMesh->getMPIGlobalCommunicator();
        MPI_Comm_rank(comm,&rank);
        MPI_Comm_size(comm,&npes);

        unsigned int activeCommSz;

        char fName[256];
        unsigned int restoreStatus=0;
        unsigned int restoreStatusGlobal=0; // 0 indicates successfully restorable.

        ot::Mesh* newMesh;
        unsigned int restoreStep[2];
        restoreStep[0]=0;
        restoreStep[1]=0;

        unsigned int restoreFileIndex=0;

        for(unsigned int cpIndex=0; cpIndex<2; cpIndex++) {

            restoreStatus=0;

            if(!rank)
            {
                sprintf(fName,"%s_step_%d.cp",EM3_CHKPT_FILE_PREFIX.c_str() , cpIndex);
                std::ifstream infile(fName);
                if(!infile) {
                    std::cout<<fName<<" file open failed "<<std::endl;
                    restoreStatus=1;
                }


                if(restoreStatus==0)
                {
                    infile>>checkPoint;
                    m_uiTinfo._m_uiTb    = checkPoint["DENDRO_TS_TIME_BEGIN"];
                    m_uiTinfo._m_uiTe    = checkPoint["DENDRO_TS_TIME_END"];
                    m_uiTinfo._m_uiT     = checkPoint["DENDRO_TS_TIME_CURRENT"];
                    m_uiTinfo._m_uiStep  = checkPoint["DENDRO_TS_STEP_CURRENT"]; 
                    m_uiTinfo._m_uiTh    = checkPoint["DENDRO_TS_TIME_STEP_SIZE"];
                    m_uiElementOrder     = checkPoint["DENDRO_TS_ELEMENT_ORDER"];
                    
                    EM3_WAVELET_TOL=checkPoint["DENDRO_TS_WAVELET_TOLERANCE"];
                    EM3_LOAD_IMB_TOL=checkPoint["DENDRO_TS_LOAD_IMB_TOLERANCE"];
                    
                    numVars=checkPoint["DENDRO_TS_NUM_VARS"];
                    activeCommSz=checkPoint["DENDRO_TS_ACTIVE_COMM_SZ"];
                    
                    restoreStep[cpIndex]=m_uiTinfo._m_uiStep;

                }
              }

            }

            if(!rank)
            {
                if(restoreStep[0]<restoreStep[1])
                    restoreFileIndex=1;
                else
                    restoreFileIndex=0;
            }

            par::Mpi_Bcast(&restoreFileIndex,1,0,comm);

            restoreStatus=0;
            octree.clear();
            if(!rank) std::cout<<"[BSSNCtx] :  Trying to restore from checkpoint index : "<<restoreFileIndex<<std::endl;
        
            if(!rank)
            {
                sprintf(fName,"%s_step_%d.cp", EM3_CHKPT_FILE_PREFIX.c_str(), restoreFileIndex);
                std::ifstream infile(fName);
                if(!infile) {
                    std::cout<<fName<<" file open failed "<<std::endl;
                    restoreStatus=1;
                }


                if(restoreStatus==0)
                {
                    infile>>checkPoint;
                    m_uiTinfo._m_uiTb    = checkPoint["DENDRO_TS_TIME_BEGIN"];
                    m_uiTinfo._m_uiTe    = checkPoint["DENDRO_TS_TIME_END"];
                    m_uiTinfo._m_uiT     = checkPoint["DENDRO_TS_TIME_CURRENT"];
                    m_uiTinfo._m_uiStep  = checkPoint["DENDRO_TS_STEP_CURRENT"]; 
                    m_uiTinfo._m_uiTh    = checkPoint["DENDRO_TS_TIME_STEP_SIZE"];
                    m_uiElementOrder     = checkPoint["DENDRO_TS_ELEMENT_ORDER"];
                    
                    EM3_WAVELET_TOL=checkPoint["DENDRO_TS_WAVELET_TOLERANCE"];
                    EM3_LOAD_IMB_TOL=checkPoint["DENDRO_TS_LOAD_IMB_TOLERANCE"];
                    
                    numVars=checkPoint["DENDRO_TS_NUM_VARS"];
                    activeCommSz=checkPoint["DENDRO_TS_ACTIVE_COMM_SZ"];
                    
                    restoreStep[restoreFileIndex]=m_uiTinfo._m_uiStep;
                
                }


            }

            par::Mpi_Allreduce(&restoreStatus,&restoreStatusGlobal,1,MPI_MAX,comm);
            if(restoreStatusGlobal==1) 
            {
                if(!rank)
                    std::cout<<"[BSSNCtx] : Restore step failed, restore file corrupted. "<<std::endl;
                MPI_Abort(comm,0);
            }


            MPI_Bcast(&m_uiTinfo,sizeof(ts::TSInfo),MPI_BYTE,0,comm);
            par::Mpi_Bcast(&EM3_WAVELET_TOL,1,0,comm);
            par::Mpi_Bcast(&EM3_LOAD_IMB_TOL,1,0,comm);

            par::Mpi_Bcast(&numVars,1,0,comm);
            par::Mpi_Bcast(&m_uiElementOrder,1,0,comm);
            par::Mpi_Bcast(&activeCommSz,1,0,comm);
            
            if(activeCommSz>npes)
            {
                if(!rank)
                    std::cout<<" [BSSNCtx] : checkpoint file written from  a larger communicator than the current global comm. (i.e. communicator shrinking not allowed in the restore step. )"<<std::endl;
                
                MPI_Abort(comm,0);
            }



            bool isActive=(rank<activeCommSz);

            MPI_Comm newComm;
            par::splitComm2way(isActive,&newComm,comm);

            if(isActive) {

                int activeRank;
                int activeNpes;

                MPI_Comm_rank(newComm, &activeRank);
                MPI_Comm_size(newComm, &activeNpes);
                assert(activeNpes == activeCommSz);

                sprintf(fName, "%s_octree_%d_%d.oct", EM3_CHKPT_FILE_PREFIX.c_str(),restoreFileIndex,activeRank);
                restoreStatus=io::checkpoint::readOctFromFile(fName, octree);
                assert(par::test::isUniqueAndSorted(octree, newComm));

            }

            par::Mpi_Allreduce(&restoreStatus,&restoreStatusGlobal,1,MPI_MAX,comm);
            if(restoreStatusGlobal==1) {

                if(!rank) std::cout<<"[BSSNCtx]: octree (*.oct) restore file is corrupted "<<std::endl;
                MPI_Abort(comm,0);
            }

            newMesh=new ot::Mesh(octree,1,m_uiElementOrder,activeCommSz,comm);
            // no need to transfer data only to resize the contex variables. 
            this->grid_transfer( newMesh , false, false, false);
            this->update_app_vars();
            
            // only reads the evolution variables. 
            if(isActive) {

                int activeRank;
                int activeNpes;

                DendroScalar** inVec=NULL;
                m_uiEVar.Get2DArray(inVec,false);

                MPI_Comm_rank(newComm, &activeRank);
                MPI_Comm_size(newComm, &activeNpes);
                assert(activeNpes == activeCommSz);

                sprintf(fName,"%s_%d_%d.var",EM3_CHKPT_FILE_PREFIX.c_str(),restoreFileIndex,activeRank);
                restoreStatus=io::checkpoint::readVecFromFile(fName,newMesh,inVec,EM3_NUM_VARS);

                delete [] inVec;
            }

            MPI_Comm_free(&newComm);
            par::Mpi_Allreduce(&restoreStatus,&restoreStatusGlobal,1,MPI_MAX,comm);
            if(restoreStatusGlobal==1) {

                if(!rank) std::cout<<"[BSSNCtx]: varible (*.var) restore file currupted "<<std::endl;
                MPI_Abort(comm,0);
            }

            std::swap(m_uiMesh,newMesh);
            delete newMesh;
            
        unsigned int localSz=m_uiMesh->getNumLocalMeshElements();
        unsigned int totalElems=0;
        par::Mpi_Allreduce(&localSz, &totalElems ,1,MPI_SUM,comm);
        
        if(!rank) std::cout<<" checkpoint at step : "<<m_uiTinfo._m_uiStep<<"active Comm. sz: "<<activeCommSz<<" restore successful: "<<" restored mesh size: "<<totalElems<<std::endl;

        m_uiIsETSSynced = false;
        return 0;

    }

    int EM3Ctx::pre_timestep(DVec sIn)
    {
        return 0;
    }

    int EM3Ctx::pre_stage(DVec  sIn)
    {
        return 0;
    }

    int EM3Ctx::post_stage(DVec sIn)
    { 
        return 0;
    }

    int EM3Ctx::post_timestep(DVec sIn)
    {
        return 0;
    }

    bool EM3Ctx::is_remesh()
    {
        
        #ifdef __PROFILE_CTX__
            m_uiCtxpt[ts::CTXPROFILE::IS_REMESH].start();
        #endif

        bool isRefine = false;

        if(EM3_ENABLE_BLOCK_ADAPTIVITY)
            return isRefine;
        
        MPI_Comm comm = m_uiMesh->getMPIGlobalCommunicator();
        this->unzip(m_uiEVar,m_uiEUnzip[0],EM3_ASYNC_COMM_K);

        DendroScalar** unzipVar;
        m_uiEUnzip[0].Get2DArray(unzipVar, false);

        unsigned int refineVarIds[EM3_NUM_REFINE_VARS];
        for(unsigned int vIndex=0; vIndex<EM3_NUM_REFINE_VARS; vIndex++)
            refineVarIds[vIndex]=EM3_REFINE_VARIABLE_INDICES[vIndex];

        double wTol=EM3_WAVELET_TOL;
        std::function<double(double,double,double,double*)> waveletTolFunc =[wTol](double x,double y, double z,double*hx) {
            return computeWTol(x,y,z,hx);
        };


        isRefine=m_uiMesh->isReMeshUnzip((const double **)unzipVar,refineVarIds,EM3_NUM_REFINE_VARS,waveletTolFunc,EM3_DENDRO_AMR_FAC); 

        delete [] unzipVar;
        return isRefine;

        #ifdef __PROFILE_CTX__
            m_uiCtxpt[ts::CTXPROFILE::IS_REMESH].stop();
        #endif

    }

    int EM3Ctx::update_app_vars() 
    {
        m_uiEVar = m_uiEvolutionVar[0];
        m_uiPVar = m_uiPrimitiveVar[0];
        m_uiCVar = m_uiConstrainedVar[0];
       
        m_uiEUnzip[0] = m_uiEvolutionUnzipVar[0];
        m_uiEUnzip[1] = m_uiEvolutionUnzipVar[1];

        m_uiCUnzip[0] = m_uiConstraintUnzipVar[0];
        m_uiCUnzip[1] = m_uiConstraintUnzipVar[1];

        return 0;
        
    }

    DVec EM3Ctx::get_evolution_vars()
    {
        return m_uiEVar;
    }
    
    DVec EM3Ctx::get_constraint_vars()
    {
        return m_uiCVar;
    }

    DVec EM3Ctx::get_primitive_vars()
    {
        return m_uiPVar;
    }

    int EM3Ctx::compute_constraints()
    {   

        if(!m_uiMesh->isActive())
            return 0;

        DendroScalar * consVar[EM3_CONSTRAINT_NUM_VARS];
        DendroScalar * evolVar[EM3_NUM_VARS];
        DendroScalar * primVar[EM3_NUM_VARS];

        m_uiEVar.Get2DVec(evolVar);
        m_uiPVar.Get2DVec(primVar);
        m_uiCVar.Get2DVec(consVar);


        DendroScalar * eunzipIn[em3::EM3_NUM_VARS];
        DendroScalar * cunzipIn[em3::EM3_CONSTRAINT_NUM_VARS];
        DendroScalar * cunzipOut[em3::EM3_CONSTRAINT_NUM_VARS]; 
    
        m_uiEUnzip[0].Get2DVec(eunzipIn);
        m_uiCUnzip[0].Get2DVec(cunzipIn);
        m_uiCUnzip[1].Get2DVec(cunzipOut);

        this->unzip(m_uiEVar,m_uiEUnzip[0],this->get_async_batch_sz());
            
        const std::vector<ot::Block> blkList=m_uiMesh->getLocalBlockList();
        unsigned int offset;
        double ptmin[3], ptmax[3];
        unsigned int sz[3];
        unsigned int bflag;
        double dx,dy,dz;
        const Point pt_min(em3::EM3_COMPD_MIN[0],em3::EM3_COMPD_MIN[1],em3::EM3_COMPD_MIN[2]);
        const Point pt_max(em3::EM3_COMPD_MAX[0],em3::EM3_COMPD_MAX[1],em3::EM3_COMPD_MAX[2]);
        const unsigned int PW=em3::EM3_PADDING_WIDTH;
        
        for(unsigned int blk=0; blk<blkList.size(); blk++)
        {
            offset=blkList[blk].getOffset();
            sz[0]=blkList[blk].getAllocationSzX();
            sz[1]=blkList[blk].getAllocationSzY();
            sz[2]=blkList[blk].getAllocationSzZ();

            bflag=blkList[blk].getBlkNodeFlag();

            dx=blkList[blk].computeDx(pt_min,pt_max);
            dy=blkList[blk].computeDy(pt_min,pt_max);
            dz=blkList[blk].computeDz(pt_min,pt_max);

            ptmin[0]=GRIDX_TO_X(blkList[blk].getBlockNode().minX())-PW*dx;
            ptmin[1]=GRIDY_TO_Y(blkList[blk].getBlockNode().minY())-PW*dy;
            ptmin[2]=GRIDZ_TO_Z(blkList[blk].getBlockNode().minZ())-PW*dz;

            ptmax[0]=GRIDX_TO_X(blkList[blk].getBlockNode().maxX())+PW*dx;
            ptmax[1]=GRIDY_TO_Y(blkList[blk].getBlockNode().maxY())+PW*dy;
            ptmax[2]=GRIDZ_TO_Z(blkList[blk].getBlockNode().maxZ())+PW*dz;

            physical_constraints(cunzipOut, (const DendroScalar **) eunzipIn, offset, ptmin, ptmax, sz, bflag);
        }

        this->zip(m_uiCUnzip[1],m_uiCVar,this->get_async_batch_sz());
        
        DendroIntL local_dof=m_uiMesh->getNumLocalMeshNodes();
        DendroIntL total_dof=0;
        par::Mpi_Reduce(&local_dof,&total_dof,1,MPI_SUM,0,m_uiMesh->getMPICommunicator()); 

        DendroScalar l2_norm[EM3_CONSTRAINT_NUM_VARS];
        DendroScalar l2_rs  [EM3_CONSTRAINT_NUM_VARS];
        DendroScalar vmin   [EM3_CONSTRAINT_NUM_VARS];
        DendroScalar vmax   [EM3_CONSTRAINT_NUM_VARS];

        for(unsigned int v=0; v < EM3_CONSTRAINT_NUM_VARS; v++)
        {
            l2_norm[v] = normL2(m_uiMesh,consVar[v],ot::VEC_TYPE::CG_NODAL,true);
            l2_rs[v]   = rsNormLp(m_uiMesh,consVar[v],2);
            vmin[v]    = vecMin(m_uiMesh,consVar[v],ot::VEC_TYPE::CG_NODAL,true); 
            vmax[v]    = vecMax(m_uiMesh,consVar[v],ot::VEC_TYPE::CG_NODAL,true); 
        }
        

        if(!m_uiMesh->getMPIRank()) {
            std::cout << "executing step: " << m_uiTinfo._m_uiStep << " dt: " << m_uiTinfo._m_uiTh << " rk_time : "<< m_uiTinfo._m_uiT << std::endl;

            for(unsigned int v=0; v <EM3_CONSTRAINT_NUM_VARS; v++)
                std::cout << "\t ||"<<EM3_CONSTRAINT_VAR_NAMES[v]<<"|| (min, max, l2,l2rs) : ("<<vmin[v]<<", "<<vmax[v]<<", "<<l2_norm[v]<<", "<<l2_rs[v]<<" ) "<<std::endl;

        }

        if(!m_uiMesh->getMPIRank())
        {

            std::ofstream fileGW;
            char fName[256];
            sprintf(fName,"%s_Constraints.dat",EM3_PROFILE_FILE_PREFIX.c_str());
            fileGW.open (fName,std::ofstream::app);
            
            // writes the header
            if(m_uiTinfo._m_uiStep==0)
            {
                fileGW<<"timestep\t"
                      <<"time\t"
                      <<"gridpts\t";
                
                for(unsigned int v=0; v < EM3_CONSTRAINT_NUM_VARS; v++)
                    fileGW<<EM3_CONSTRAINT_VAR_NAMES[v]<<"_l2\t"<<EM3_CONSTRAINT_VAR_NAMES[v]<<"_l2rs\t";
                      
                fileGW<<std::endl;
            }
                

            fileGW<<m_uiTinfo._m_uiStep<<"\t"
                  <<m_uiTinfo._m_uiT   <<"\t"
                  <<total_dof          <<"\t";

            for(unsigned int v=0; v < EM3_CONSTRAINT_NUM_VARS; v++)
                fileGW<<l2_norm[v]<<"\t"<<l2_rs[v]<<"\t";
            
            fileGW<<std::endl;
            fileGW.close();

        }



        
        return 0;
    }

    int EM3Ctx::compute_primitives()
    {

        if(!m_uiMesh->isActive())
            return 0; 
        

        DendroScalar * consVar[EM3_CONSTRAINT_NUM_VARS];
        DendroScalar * evolVar[EM3_NUM_VARS];
        DendroScalar * primVar[EM3_NUM_VARS];

        m_uiEVar.Get2DVec(evolVar);
        m_uiPVar.Get2DVec(primVar);
        m_uiCVar.Get2DVec(consVar);

        // initialize diff begin.
        unsigned int nodeLookUp_CG;
        unsigned int nodeLookUp_DG;
        double x,y,z,len;
        const ot::TreeNode * pNodes=&(*(m_uiMesh->getAllElements().begin()));
        unsigned int ownerID,ii_x,jj_y,kk_z;
        unsigned int eleOrder=m_uiMesh->getElementOrder();
        const unsigned int * e2n_cg=&(*(m_uiMesh->getE2NMapping().begin()));
        const unsigned int * e2n_dg=&(*(m_uiMesh->getE2NMapping_DG().begin()));
        const unsigned int nPe=m_uiMesh->getNumNodesPerElement();
        const unsigned int nodeLocalBegin=m_uiMesh->getNodeLocalBegin();
        const unsigned int nodeLocalEnd=m_uiMesh->getNodeLocalEnd();

        double var[em3::EM3_NUM_VARS];
        for(unsigned int elem=m_uiMesh->getElementLocalBegin();elem<m_uiMesh->getElementLocalEnd();elem++)
        {

            for(unsigned int k=0;k<(eleOrder+1);k++)
                for(unsigned int j=0;j<(eleOrder+1);j++ )
                    for(unsigned int i=0;i<(eleOrder+1);i++)
                    {
                        nodeLookUp_CG=e2n_cg[elem*nPe+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i];
                        if(nodeLookUp_CG>=nodeLocalBegin && nodeLookUp_CG<nodeLocalEnd)
                        {
                            nodeLookUp_DG=e2n_dg[elem*nPe+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i];
                            m_uiMesh->dg2eijk(nodeLookUp_DG,ownerID,ii_x,jj_y,kk_z);
                            len= (double) (1u<<(m_uiMaxDepth-pNodes[ownerID].getLevel()));
                            x=pNodes[ownerID].getX()+ ii_x*(len/(eleOrder));
                            y=pNodes[ownerID].getY()+ jj_y*(len/(eleOrder));
                            z=pNodes[ownerID].getZ()+ kk_z*(len/(eleOrder));
                            
                            em3::analyticalSol((double)x,(double)y,(double)z,m_uiTinfo._m_uiT,var);
                            for(unsigned int v=0; v < em3::EM3_NUM_VARS; v++)
                                primVar[v][nodeLookUp_CG] = evolVar[v][nodeLookUp_CG] -var[v];

                        }

                    }
        }

        m_uiMesh->readFromGhostBegin(m_uiPVar.GetVecArray(),m_uiPVar.GetDof());
        DendroIntL local_dof=m_uiMesh->getNumLocalMeshNodes();
        DendroIntL total_dof=0;
        par::Mpi_Reduce(&local_dof,&total_dof,1,MPI_SUM,0,m_uiMesh->getMPICommunicator()); 

        double l2_norm[em3::EM3_NUM_VARS];
        for(unsigned int v=0; v < em3::EM3_NUM_VARS; v++)
            l2_norm[v] = normL2(m_uiMesh,primVar[v],ot::VEC_TYPE::CG_NODAL,true);

        if(!m_uiMesh->getMPIRank())
        {

            std::ofstream fileGW;
            char fName[256];
            sprintf(fName,"%s_Primitives.dat",em3::EM3_PROFILE_FILE_PREFIX.c_str());
            fileGW.open (fName,std::ofstream::app);
            // writes the header
            if(m_uiTinfo._m_uiStep==0)
            {  
                fileGW  <<"timestep\t"
                        <<"time\t"
                        <<"gridpts\t";

                for(unsigned int v=0; v < EM3_NUM_VARS; v++)
                    fileGW<<"diff_"<<EM3_VAR_NAMES[v]<<"_l2\t";

                fileGW<<std::endl;

            }

            fileGW<<m_uiTinfo._m_uiStep<<"\t"
                  <<m_uiTinfo._m_uiT   <<"\t"
                  <<total_dof          <<"\t";

                  for(unsigned int v=0; v < EM3_NUM_VARS; v++)
                    fileGW<<l2_norm[v]<<"\t";

                  fileGW<<std::endl;

            fileGW.close();
            

        }
        m_uiMesh->readFromGhostEnd(m_uiPVar.GetVecArray(),m_uiPVar.GetDof());
        return 0;

    }

    int EM3Ctx::terminal_output()
    {
        if(m_uiMesh->isActive())
        {
            
            const unsigned int currentStep = m_uiTinfo._m_uiStep;

            if((m_uiMesh->isActive()) && (currentStep%em3::EM3_TIME_STEP_OUTPUT_FREQ)==0)
            {
                for(unsigned int v=0; v < m_uiEVar.GetDof(); v++)
                {
                    DendroScalar min=0, max=0;
                    m_uiEVar.VecMinMax(m_uiMesh,min,max,v);
                    if(!(m_uiMesh->getMPIRank()))
                    std::cout<<"[EM3Ctx]:  "<<EM3_VAR_NAMES[v]<<" (min,max) : \t ( "<<min<<", "<<max<<" ) "<<std::endl;
                }

                // compare with the analytical solution. 
                DendroScalar * consVar[EM3_CONSTRAINT_NUM_VARS];
                DendroScalar * evolVar[EM3_NUM_VARS];
                DendroScalar * primVar[EM3_NUM_VARS];

                m_uiEVar.Get2DVec(evolVar);
                m_uiPVar.Get2DVec(primVar);
                m_uiCVar.Get2DVec(consVar);
                
                this->compute_primitives();
                
                DendroScalar l2_norm[EM3_NUM_VARS];
                DendroScalar l2_rs  [EM3_NUM_VARS];
                DendroScalar vmin   [EM3_NUM_VARS];
                DendroScalar vmax   [EM3_NUM_VARS];

                for(unsigned int v=0; v < EM3_NUM_VARS; v++)
                {
                    l2_norm[v] = normL2(m_uiMesh,primVar[v],ot::VEC_TYPE::CG_NODAL,true);
                    l2_rs[v]   = rsNormLp(m_uiMesh,primVar[v],2);
                    vmin[v]    = vecMin(m_uiMesh,primVar[v],ot::VEC_TYPE::CG_NODAL,true); 
                    vmax[v]    = vecMax(m_uiMesh,primVar[v],ot::VEC_TYPE::CG_NODAL,true); 
                }
                

                if(!m_uiMesh->getMPIRank()) {
                    std::cout << "executing step: " << m_uiTinfo._m_uiStep << " dt: " << m_uiTinfo._m_uiTh << " rk_time : "<< m_uiTinfo._m_uiT << std::endl;

                    for(unsigned int v=0; v <EM3_NUM_VARS; v++)
                        std::cout << "\t ||diff_["<<EM3_VAR_NAMES[v]<<"]|| (min, max, l2,l2rs) : ("<<vmin[v]<<", "<<vmax[v]<<", "<<l2_norm[v]<<", "<<l2_rs[v]<<" ) "<<std::endl;

                }


            }

        }


        return 0;
    }

    unsigned int EM3Ctx::getBlkTimestepFac(unsigned int blev, unsigned int lmin, unsigned int lmax)
    {
        const unsigned int ldiff=0;
        if((lmax-blev) <=ldiff)
            return 1;
        else
        {
            return 1u<<(lmax-blev-ldiff);
        }

    }



} // end of namespace em3