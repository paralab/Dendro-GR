//
// Created by milinda on 8/10/18.
//
#include "rhs_cuda.cuh"


namespace cuda
{



    void computeRHS(double **unzipVarsRHS, const double **uZipVars,const ot::Block* dendroBlockList,unsigned int numBlocks,const cuda::BSSNComputeParams* bssnPars,dim3 blockDim,const Point & pt_min, const Point & pt_max,unsigned int numStreams,unsigned int device)
    {
        cuda::profile::t_overall.start();


        if(numStreams==0)
        {
            std::cout<<"[Error]: "<<__func__<<" numStreams "<<numStreams<<" should at least be 1 (synchronous transfer) "<<std::endl;
            return ;
        }

        // initialize the input data
        unsigned int offset,bflag;
        unsigned int sz[3];
        double dx,dy,dz;
        double ptmin[3];//={-1.0,-1.0,-1.0};
        double ptmax[3];
        double hx[3];


        unsigned int c_in=0;
        unsigned int c_ex=0;

        for(unsigned int blk=0;blk<numBlocks;blk++)
        {
            if(dendroBlockList[blk].getBlkNodeFlag()==0)
                c_in++;
            else
                c_ex++;

        }

        const unsigned int NUM_GPU_DENDRO_BLOCKS=c_in;
        const unsigned int NUM_CPU_DENDRO_BLOCKS=c_ex;
        const unsigned int BSSN_NUM_VARS=24;
        const unsigned int BSSN_CONSTRAINT_NUM_VARS=6;
        const unsigned int UNZIP_DOF_SZ=dendroBlockList[numBlocks-1].getOffset()+ dendroBlockList[numBlocks-1].getAlignedBlockSz();



        //get GPU information.
        // assumes the if there are multiple gpus per node all have the same specification.
        cudaSetDevice(device);
        cuda::__CUDA_DEVICE_PROPERTIES=getGPUDeviceInfo(device);
        // device properties for the host
        cudaDeviceProp deviceProp;
        cudaGetDeviceProperties(&deviceProp,device);

        const double GPU_BLOCK_SHARED_MEM_UTIL=0.8;


        ot::Block* blkListCPU=new ot::Block[NUM_CPU_DENDRO_BLOCKS];
        ot::Block* blkListGPU=new ot::Block[NUM_GPU_DENDRO_BLOCKS];
        cuda::_Block* cudaBlkList=new cuda::_Block[NUM_GPU_DENDRO_BLOCKS];

        c_in=0;
        c_ex=0;

        for(unsigned int blk=0;blk<numBlocks;blk++)
        {

            bflag=dendroBlockList[blk].getBlkNodeFlag();
            if(bflag==0)
            {
                offset=dendroBlockList[blk].getOffset();
                sz[0]=dendroBlockList[blk].getAllocationSzX();
                sz[1]=dendroBlockList[blk].getAllocationSzY();
                sz[2]=dendroBlockList[blk].getAllocationSzZ();



                dx = dendroBlockList[blk].computeDx(pt_min,pt_max);
                dy = dendroBlockList[blk].computeDy(pt_min,pt_max);
                dz = dendroBlockList[blk].computeDz(pt_min,pt_max);

                hx[0]=dx;
                hx[1]=dy;
                hx[2]=dz;


                ptmin[0]=GRIDX_TO_X(dendroBlockList[blk].getBlockNode().minX())-3*dx;
                ptmin[1]=GRIDY_TO_Y(dendroBlockList[blk].getBlockNode().minY())-3*dy;
                ptmin[2]=GRIDZ_TO_Z(dendroBlockList[blk].getBlockNode().minZ())-3*dz;


                ptmax[0]=GRIDX_TO_X(dendroBlockList[blk].getBlockNode().maxX())+3*dx;
                ptmax[1]=GRIDY_TO_Y(dendroBlockList[blk].getBlockNode().maxY())+3*dy;
                ptmax[2]=GRIDZ_TO_Z(dendroBlockList[blk].getBlockNode().maxZ())+3*dz;

                cudaBlkList[c_in]=cuda::_Block((const double *)ptmin,(const double *)ptmax,offset,bflag,(const unsigned int*)sz, (const double *)hx);
                blkListGPU[c_in]=dendroBlockList[blk];
                c_in++;

            }else {
                blkListCPU[c_ex]=dendroBlockList[blk];
                c_ex++;
            }

        }


        unsigned int maxBlkSz=0;
        for(unsigned int blk=0;blk<NUM_GPU_DENDRO_BLOCKS;blk++)
        {
            if(maxBlkSz<cudaBlkList[blk].getAlignedBlockSz())
                maxBlkSz=cudaBlkList[blk].getAlignedBlockSz();
        }

        const unsigned int derivSz=(maxBlkSz);
        cuda::__DENDRO_BLK_MAX_SZ=cuda::copyValueToDevice(&derivSz);
        const unsigned int numSM=deviceProp.multiProcessorCount;

        //std::cout<<"deriv alloc begin"<<std::endl;

        cuda::profile::t_cudaMalloc_derivs.start();

        cuda::MemoryDerivs derivWorkSpace;
        derivWorkSpace.allocateDerivMemory(maxBlkSz,numSM,numStreams);
        CUDA_CHECK_ERROR();

        cuda::__BSSN_DERIV_WORKSPACE=cuda::copyValueToDevice(&derivWorkSpace);
        CUDA_CHECK_ERROR();

        cuda::profile::t_cudaMalloc_derivs.stop();


        if(numStreams==1)
        {// sync case
            // computes 1D grid block
            dim3 gridDim;
            unsigned int *blockMap=NULL;
            cuda::computeDendroBlockToGPUMap(blkListGPU,NUM_GPU_DENDRO_BLOCKS,blockMap,gridDim);

            cuda::profile::t_H2D_Comm.start();
            const unsigned int NUM_GPU_BLOCKS=((gridDim.x)*(gridDim.y)*(gridDim.z));




            //send blocks to the gpu
            cuda::__DENDRO_BLOCK_LIST=cuda::copyArrayToDevice(cudaBlkList,NUM_GPU_DENDRO_BLOCKS);
            cuda::__DENDRO_NUM_BLOCKS=cuda::copyValueToDevice(&NUM_GPU_DENDRO_BLOCKS);

            cuda::__NUM_GPU_BLOCKS=cuda::copyValueToDevice(&NUM_GPU_BLOCKS);
            cuda::__GPU_BLOCK_MAP=cuda::copyArrayToDevice(blockMap,2*NUM_GPU_BLOCKS);

            cuda::__BSSN_NUM_VARS=cuda::copyValueToDevice(&BSSN_NUM_VARS);
            cuda::__BSSN_CONSTRAINT_NUM_VARS=cuda::copyValueToDevice(&BSSN_CONSTRAINT_NUM_VARS);

            cuda::__GPU_BLOCK_SHARED_MEM_UTIL=cuda::copyValueToDevice(&GPU_BLOCK_SHARED_MEM_UTIL);

            //allocate memory for unzip vectors
            cuda::__UNZIP_INPUT=cuda::alloc2DCudaArray<double>(uZipVars,BSSN_NUM_VARS,UNZIP_DOF_SZ);
            cuda::__UNZIP_OUTPUT=cuda::alloc2DCudaArray<double>(BSSN_NUM_VARS,UNZIP_DOF_SZ);

            cuda::__BSSN_COMPUTE_PARMS=cuda::copyValueToDevice(&(*bssnPars));


            cuda::profile::t_H2D_Comm.stop();


            cuda::profile::t_rhs_total.start();


            cuda::profile::t_rhs_gpu.start();

            cuda::__computeBSSNRHS<<<gridDim,blockDim>>>(cuda::__UNZIP_OUTPUT,(const double**)cuda::__UNZIP_INPUT,cuda::__BSSN_DERIV_WORKSPACE,cuda::__DENDRO_BLOCK_LIST,cuda::__GPU_BLOCK_MAP,cuda::__BSSN_COMPUTE_PARMS,cuda::__CUDA_DEVICE_PROPERTIES,0);

            cudaDeviceSynchronize();
            CUDA_CHECK_ERROR();

            cuda::profile::t_rhs_gpu.stop();


            cuda::profile::t_D2H_Comm.start();
            cuda::copy2DArrayToHost(unzipVarsRHS,(const double**)cuda::__UNZIP_OUTPUT,BSSN_NUM_VARS,UNZIP_DOF_SZ);
            cuda::profile::t_D2H_Comm.stop();


            cuda::profile::t_rhs_cpu.start();

            for(unsigned int blk=0;blk<NUM_CPU_DENDRO_BLOCKS;blk++)
            {
                offset=blkListCPU[blk].getOffset();
                sz[0]=blkListCPU[blk].getAllocationSzX();
                sz[1]=blkListCPU[blk].getAllocationSzY();
                sz[2]=blkListCPU[blk].getAllocationSzZ();

                bflag=blkListCPU[blk].getBlkNodeFlag();

                dx=blkListCPU[blk].computeDx(pt_min,pt_max);
                dy=blkListCPU[blk].computeDy(pt_min,pt_max);
                dz=blkListCPU[blk].computeDz(pt_min,pt_max);

                ptmin[0]=GRIDX_TO_X(blkListCPU[blk].getBlockNode().minX())-3*dx;
                ptmin[1]=GRIDY_TO_Y(blkListCPU[blk].getBlockNode().minY())-3*dy;
                ptmin[2]=GRIDZ_TO_Z(blkListCPU[blk].getBlockNode().minZ())-3*dz;

                ptmax[0]=GRIDX_TO_X(blkListCPU[blk].getBlockNode().maxX())+3*dx;
                ptmax[1]=GRIDY_TO_Y(blkListCPU[blk].getBlockNode().maxY())+3*dy;
                ptmax[2]=GRIDZ_TO_Z(blkListCPU[blk].getBlockNode().maxZ())+3*dz;

                bssnrhs(unzipVarsRHS, (const double **)uZipVars, offset, ptmin, ptmax, sz, bflag);

            }

            cuda::profile::t_rhs_cpu.stop();


            cuda::profile::t_rhs_total.stop();

            delete [] blockMap;

        }else {

            cudaStream_t streams[numStreams];
            // send counts and offset for asynctransfer, all the counts and offsets are based on the how we transfer the data blocks to the gpu.
            unsigned int sendBlockCount[numStreams];
            unsigned int sendBlockOffset[numStreams];

            unsigned int sendUnzipCount[numStreams];
            unsigned int sendUnzipOffset[numStreams];

            unsigned int gridBlockCount[numStreams];
            unsigned int gridBlockOffset[numStreams];

            for(unsigned int i=0;i<numStreams;i++)
            {
                cudaStreamCreate(&streams[i]);
                CUDA_CHECK_ERROR();
            }


            for(unsigned int i=0;i<numStreams;i++)
            {
                sendBlockCount[i]=(((i+1)*NUM_GPU_DENDRO_BLOCKS)/numStreams)-(((i)*NUM_GPU_DENDRO_BLOCKS)/numStreams);
                sendUnzipCount[i]=0;
                for(unsigned int blk=(((i)*NUM_GPU_DENDRO_BLOCKS)/numStreams);blk<(((i+1)*NUM_GPU_DENDRO_BLOCKS)/numStreams);blk++)
                    sendUnzipCount[i]+=blkListGPU[blk].getAlignedBlockSz();

            }


            sendBlockOffset[0]=0;
            sendUnzipOffset[0]=0;
            for(unsigned int i=1;i<numStreams;i++)
            {
                sendBlockOffset[i]=sendBlockOffset[i-1]+sendBlockCount[i-1];
                sendUnzipOffset[i]=sendUnzipOffset[i-1]+sendUnzipCount[i-1];
            }


            unsigned int ** blockMap=new unsigned int*[numStreams];
            dim3 gridDim[numStreams];

            for(unsigned int i=0;i<numStreams;i++)
            {
                cuda::computeDendroBlockToGPUMap(blkListGPU+sendBlockOffset[i],sendBlockCount[i],blockMap[i],gridDim[i]);
                gridBlockCount[i]=(gridDim[i].x*gridDim[i].y*gridDim[i].z);

                // modify for map for the global blk list index.
                for(unsigned int gblk=0;gblk<gridBlockCount[i];gblk++)
                {
                    blockMap[i][2*gblk]+=sendBlockOffset[i];
                    blockMap[i][2*gblk+1]+=sendBlockOffset[i];

                }

            }

            gridBlockOffset[0]=0;
            for(unsigned int i=1;i<numStreams;i++)
                gridBlockOffset[i]=gridBlockOffset[i-1]+gridBlockCount[i-1];


            const unsigned int NUM_GPU_BLOCKS=gridBlockOffset[numStreams-1]+gridBlockCount[numStreams-1];
            double ** unZipInputHost=NULL;
            double ** unZipOutputHost=NULL;

            cuda::profile::t_H2D_Comm.start();

            cuda::__DENDRO_BLOCK_LIST=cuda::alloc1DCudaArray<_Block>(NUM_GPU_DENDRO_BLOCKS);
            cuda::__GPU_BLOCK_MAP=cuda::alloc1DCudaArray<unsigned int>(2*NUM_GPU_BLOCKS);
            cuda::__UNZIP_INPUT=cuda::alloc2DCudaArray<double>(unZipInputHost,BSSN_NUM_VARS,UNZIP_DOF_SZ);
            cuda::__UNZIP_OUTPUT=cuda::alloc2DCudaArray<double>(unZipOutputHost,BSSN_NUM_VARS,UNZIP_DOF_SZ);

            // all the blocking small sends

            cuda::__DENDRO_NUM_BLOCKS=cuda::copyValueToDevice(&NUM_GPU_DENDRO_BLOCKS);
            cuda::__NUM_GPU_BLOCKS=cuda::copyValueToDevice(&NUM_GPU_BLOCKS);
            cuda::__BSSN_NUM_VARS=cuda::copyValueToDevice(&BSSN_NUM_VARS);
            cuda::__BSSN_CONSTRAINT_NUM_VARS=cuda::copyValueToDevice(&BSSN_CONSTRAINT_NUM_VARS);
            cuda::__BSSN_COMPUTE_PARMS=cuda::copyValueToDevice(&(*bssnPars));

            cuda::profile::t_H2D_Comm.stop();


            cuda::profile::t_rhs_total.start();


            for(unsigned int i=0;i<numStreams;i++)
            {
                for (unsigned int var = 0; var < BSSN_NUM_VARS; var++)
                    cuda::copyArrayToDeviceAsync(uZipVars[var] + sendUnzipOffset[i],unZipInputHost[var] + sendUnzipOffset[i], sendUnzipCount[i], streams[i]);


                cuda::copyArrayToDeviceAsync(cudaBlkList + sendBlockOffset[i], cuda::__DENDRO_BLOCK_LIST + sendBlockOffset[i], sendBlockCount[i], streams[i]);
                cuda::copyArrayToDeviceAsync(blockMap[i], cuda::__GPU_BLOCK_MAP + 2*gridBlockOffset[i], 2*gridBlockCount[i],streams[i]);
            }


            for(unsigned int i=0;i<numStreams;i++)
            {

                cuda::__computeBSSNRHS<<< gridDim[i], blockDim, 0, streams[i] >>>(cuda::__UNZIP_OUTPUT,(const double**)cuda::__UNZIP_INPUT,cuda::__BSSN_DERIV_WORKSPACE,cuda::__DENDRO_BLOCK_LIST,cuda::__GPU_BLOCK_MAP + 2*gridBlockOffset[i],cuda::__BSSN_COMPUTE_PARMS,cuda::__CUDA_DEVICE_PROPERTIES,i);
                CUDA_CHECK_ERROR();

            }

            for(unsigned int i=0;i<numStreams;i++)
            {
                //D2H
                for (unsigned int var = 0; var < BSSN_NUM_VARS; var++)
                    cuda::copyArrayToHostAsync(unzipVarsRHS[var]+sendUnzipOffset[i],unZipOutputHost[var]+sendUnzipOffset[i],sendUnzipCount[i],streams[i]);
            }

            cuda::profile::t_rhs_cpu.start();

            // process cpu blocks before sync
            double ** unZipRHS_CPU=new double*[BSSN_NUM_VARS];
            for(unsigned int var=0;var<BSSN_NUM_VARS;var++)
                unZipRHS_CPU[var]=new double[UNZIP_DOF_SZ];

            for(unsigned int blk=0;blk<NUM_CPU_DENDRO_BLOCKS;blk++)
            {
                offset=blkListCPU[blk].getOffset();
                sz[0]=blkListCPU[blk].getAllocationSzX();
                sz[1]=blkListCPU[blk].getAllocationSzY();
                sz[2]=blkListCPU[blk].getAllocationSzZ();

                bflag=blkListCPU[blk].getBlkNodeFlag();

                dx=blkListCPU[blk].computeDx(pt_min,pt_max);
                dy=blkListCPU[blk].computeDy(pt_min,pt_max);
                dz=blkListCPU[blk].computeDz(pt_min,pt_max);

                ptmin[0]=GRIDX_TO_X(blkListCPU[blk].getBlockNode().minX())-3*dx;
                ptmin[1]=GRIDY_TO_Y(blkListCPU[blk].getBlockNode().minY())-3*dy;
                ptmin[2]=GRIDZ_TO_Z(blkListCPU[blk].getBlockNode().minZ())-3*dz;

                ptmax[0]=GRIDX_TO_X(blkListCPU[blk].getBlockNode().maxX())+3*dx;
                ptmax[1]=GRIDY_TO_Y(blkListCPU[blk].getBlockNode().maxY())+3*dy;
                ptmax[2]=GRIDZ_TO_Z(blkListCPU[blk].getBlockNode().maxZ())+3*dz;

                bssnrhs(unZipRHS_CPU, (const double **)uZipVars, offset, ptmin, ptmax, sz, bflag);

            }

            cuda::profile::t_rhs_cpu.stop();

            cudaDeviceSynchronize();
            CUDA_CHECK_ERROR();


            // merge cpu and gpu results.

            for(unsigned int blk=0;blk<NUM_CPU_DENDRO_BLOCKS;blk++)
            {
                offset = blkListCPU[blk].getOffset();
                sz[0] = blkListCPU[blk].getAllocationSzX();
                sz[1] = blkListCPU[blk].getAllocationSzY();
                sz[2] = blkListCPU[blk].getAllocationSzZ();

                for(unsigned int var=0;var<BSSN_NUM_VARS;var++)
                    for(unsigned int node=0;node<(sz[0]*sz[1]*sz[2]);node++)
                        unzipVarsRHS[var][offset+node]=unZipRHS_CPU[var][offset+node];


            }

            cuda::profile::t_rhs_total.stop();

            for(unsigned int var=0;var<BSSN_NUM_VARS;var++)
                delete [] unZipRHS_CPU[var];

            delete [] unZipRHS_CPU;


            for(unsigned int i=0;i<numStreams;i++)
                delete [] blockMap[i];

            delete [] blockMap;
        }

        cuda::profile::t_cudaMalloc_derivs.start();
        derivWorkSpace.deallocateDerivMemory();
        CUDA_CHECK_ERROR();
        cuda::profile::t_cudaMalloc_derivs.stop();

        cudaFree(cuda::__CUDA_DEVICE_PROPERTIES);
        cudaFree(cuda::__DENDRO_BLOCK_LIST);
        cudaFree(cuda::__DENDRO_NUM_BLOCKS);
        cudaFree(cuda::__NUM_GPU_BLOCKS);
        cudaFree(cuda::__GPU_BLOCK_MAP);

        cudaFree(cuda::__BSSN_NUM_VARS);
        cudaFree(cuda::__BSSN_CONSTRAINT_NUM_VARS);
        cudaFree(cuda::__GPU_BLOCK_SHARED_MEM_UTIL);

        cudaFree(cuda::__DENDRO_BLK_MAX_SZ);
        cudaFree(cuda::__BSSN_DERIV_WORKSPACE);
        cudaFree(cuda::__BSSN_COMPUTE_PARMS);

        cuda::dealloc2DCudaArray(cuda::__UNZIP_INPUT,BSSN_NUM_VARS);
        cuda::dealloc2DCudaArray(cuda::__UNZIP_OUTPUT,BSSN_NUM_VARS);



        delete [] blkListGPU;
        delete [] blkListCPU;
        delete [] cudaBlkList;



        cuda::profile::t_overall.stop();

    }




}
