//
// Created by milinda on 9/18/18.
//


#include "gpuTest.cuh"

namespace cuda
{

    __global__ void __compute1D_derivs(double** __unzipOutVar, const double**__unzipInVar, MemoryDerivs* __derivWorkspace, const cuda::_Block* __dendroBlkList, const BSSNComputeParams* __bssnPar, const cudaDeviceProp*__deviceProperties)
    {
        const _Block dblock=__dendroBlkList[blockIdx.x];
        const unsigned int NUM_SM_UNITS=__deviceProperties->multiProcessorCount;
        const unsigned int SM_ID=blockIdx.x%NUM_SM_UNITS;
        const unsigned int offset=dblock.getOffset();
        const unsigned int *sz=dblock.getSz();
        const double* hx=dblock.getDx();
        const double dx=hx[0];
        const double dy=hx[1];
        const double dz=hx[2];
        const double* ptmin=dblock.getPtMin();
        const double* ptmax=dblock.getPtMax();
        const unsigned int bflag=dblock.getBFlag();

        const unsigned int BSSN_NUM_VARS=24;

        const unsigned int blkSz=sz[0]*sz[1]*sz[2];

        const unsigned int tile_sz[3]={729,1,1};
        const unsigned int TILE_SZ=tile_sz[0]*tile_sz[1]*tile_sz[2];

        //allocate memory for shared deriv variables.
        __shared__ double grad_0[729];
        const unsigned int Lb = 1;// load begin bound
        const unsigned int Le = sz[0]-1;// load end bound
        const unsigned int BLK_INTERATIONS = ((Le-Lb)<tile_sz[0])? 1: ((int)ceil((double)(Le-Lb-tile_sz[0])/(tile_sz[0]-2*2)))+1;

        unsigned int ijk_lm[2];

        //allocate memory for shared unzip input.
        __shared__ double unzipVarInShared[729];
        for(unsigned int iter=0;iter<BLK_INTERATIONS;iter++)
        {
            ijk_lm[0]=max(1,(int)(1 + tile_sz[0]*iter -2*iter*2));
            ijk_lm[1]=min(ijk_lm[0]+tile_sz[0],sz[0]-1);

            for(unsigned int var=0;var<BSSN_NUM_VARS;var++)
            {
                cuda::__loadGlobalToShared1D<double>(&__unzipInVar[var][offset],(double *) unzipVarInShared,(const unsigned int *) ijk_lm,(const unsigned int *) sz,(const unsigned int *) tile_sz);
                __syncthreads();

                //_RSWS_deriv42_x_1D((double *) grad_0,(const double *) unzipVarInShared,dx, (const unsigned int *) ijk_lm , (const unsigned int *) sz , (const unsigned int *) tile_sz, 2, bflag);
//                _RSWS_deriv42_x_1D((double *) grad_0,(const double *) unzipVarInShared,dx, (const unsigned int *) ijk_lm , (const unsigned int *) sz , (const unsigned int *) tile_sz, 2, bflag);
//                _RSWS_deriv42_x_1D((double *) grad_0,(const double *) unzipVarInShared,dx, (const unsigned int *) ijk_lm , (const unsigned int *) sz , (const unsigned int *) tile_sz, 2, bflag);
//                _RSWS_deriv42_x_1D((double *) grad_0,(const double *) unzipVarInShared,dx, (const unsigned int *) ijk_lm , (const unsigned int *) sz , (const unsigned int *) tile_sz, 2, bflag);
//                _RSWS_deriv42_x_1D((double *) grad_0,(const double *) unzipVarInShared,dx, (const unsigned int *) ijk_lm , (const unsigned int *) sz , (const unsigned int *) tile_sz, 2, bflag);
//                _RSWS_deriv42_x_1D((double *) grad_0,(const double *) unzipVarInShared,dx, (const unsigned int *) ijk_lm , (const unsigned int *) sz , (const unsigned int *) tile_sz, 2, bflag);
//                _RSWS_deriv42_x_1D((double *) grad_0,(const double *) unzipVarInShared,dx, (const unsigned int *) ijk_lm , (const unsigned int *) sz , (const unsigned int *) tile_sz, 2, bflag);
//                _RSWS_deriv42_x_1D((double *) grad_0,(const double *) unzipVarInShared,dx, (const unsigned int *) ijk_lm , (const unsigned int *) sz , (const unsigned int *) tile_sz, 2, bflag);
//                _RSWS_deriv42_x_1D((double *) grad_0,(const double *) unzipVarInShared,dx, (const unsigned int *) ijk_lm , (const unsigned int *) sz , (const unsigned int *) tile_sz, 2, bflag);
//                _RSWS_deriv42_x_1D((double *) grad_0,(const double *) unzipVarInShared,dx, (const unsigned int *) ijk_lm , (const unsigned int *) sz , (const unsigned int *) tile_sz, 2, bflag);

                cuda::__storeSharedToGlobal1D<double>((double *) grad_0,&__unzipOutVar[var][offset],(const unsigned int *) ijk_lm,(const unsigned int *) sz,(const unsigned int *) tile_sz);
                __syncthreads();
            }


        }



    }

} // namespace cuda

