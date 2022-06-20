//
// Created by milinda on 8/9/18.
//

/**
 * @brief Contians cuda derivs for bssn computation
 *
 * */

#include "derivs.cuh"
#include "dendro.h"

namespace cuda
{




    __device__ void _RSWS_deriv42_x_1D(double * const  Dxu, const double * const  u, const double dx, const unsigned int* ijk_lm, const unsigned int * sz, const unsigned int* tile_sz, unsigned int pw ,unsigned bflag)
    {

        const unsigned int i_b=ijk_lm[0]+pw;
        const unsigned int i_e=ijk_lm[1]-pw;


        const double idx = 1.0 / dx;
        const double idx_by_2 = 0.5 * idx;
        const double idx_by_12 = idx / 12.0;

        const int ib = 3;
        const int jb = 1;
        const int kb = 1;
        const int ie = sz[0] - 3;
        const int je = sz[1] - 1;
        const int ke = sz[2] - 1;


        for(unsigned int i=(i_b-ijk_lm[0])+threadIdx.x;i<(i_e-ijk_lm[0]);i+=blockDim.x)
        {
            Dxu[i] = (u[i - 2] - 8.0 * u[i - 1] + 8.0 * u[i + 1] - u[i + 2]) * idx_by_12;
        }


    }


/*----------------------------------------------------------------------;
 *
 * compute first derivative in x direction
 *
 *----------------------------------------------------------------------*/

    __device__ void _RSWS_deriv42_x(double * const  Dxu, const double * const  u, const double dx, const unsigned int* ijk_lm, const unsigned int * sz, const unsigned int* tile_sz,unsigned int pw, unsigned bflag) {

            const unsigned int i_b=ijk_lm[2*0+0]+pw;
            const unsigned int i_e=ijk_lm[2*0+1]-pw;

            const unsigned int j_b=max((int)ijk_lm[2*1+0],(int)1);
            const unsigned int j_e=min((int)ijk_lm[2*1+1],sz[1]-1);

            const unsigned int k_b=max((int)ijk_lm[2*2+0],int(1));
            const unsigned int k_e=min((int)ijk_lm[2*2+1],sz[2]-1);

            const unsigned int k=(k_b-ijk_lm[4])+threadIdx.y+threadIdx.z*blockDim.y;
            const unsigned int j=(j_b-ijk_lm[2])+threadIdx.x;

            if((j>(j_e-ijk_lm[2])) || (k>(k_e-ijk_lm[4])) ) return;

            const double idx = 1.0 / dx;
            const double idx_by_2 = 0.5 * idx;
            const double idx_by_12 = idx / 12.0;


            //printf("dx threadid (%d,%d,%d) loop begin: (%d,%d,%d) loop end: (%d,%d,%d)  tile begin: (%d,%d,%d) tile end: (%d,%d,%d) \n", threadIdx.x,threadIdx.y,threadIdx.z,ix_b,jy_b,kz_b,ix_e,jy_e,kz_e,ijk_lm[0],ijk_lm[2],ijk_lm[4],ijk_lm[1],ijk_lm[3],ijk_lm[5]);
            unsigned int pp = IDX_L((i_b-ijk_lm[0]), j, k);

            for(unsigned int i=(i_b-ijk_lm[0]);i<=(i_e-ijk_lm[0]);++i, ++pp)
                  Dxu[pp] = (u[pp - 2] - 8.0 * u[pp - 1] + 8.0 * u[pp + 1] - u[pp + 2]) * idx_by_12;



           /* if ((bflag & (1u << OCT_DIR_LEFT)) && (i_b==ib && (ix_b+ijk_lm[0])==ib)  ) {

                    for(unsigned int k=kz_b;k<kz_e;k++)
                        for(unsigned int j=jy_b;j<jy_e;j++)
                        {
                            Dxu[IDX_L(ix_b, j, k)] = (-3.0 * u[IDX_L(ix_b, j, k)]
                                         + 4.0 * u[IDX_L(ix_b+1, j, k)]
                                         - u[IDX_L(ix_b+2, j, k)]
                                        ) * idx_by_2;

                            Dxu[IDX_L(ix_b+1, j, k)] = (-u[IDX_L(ix_b, j, k)]
                                         + u[IDX_L(ix_b+2, j, k)]
                                        ) * idx_by_2;

                        }


                }

            if ((bflag & (1u << OCT_DIR_RIGHT)) && (i_e==ie && (ix_e+ijk_lm[0])==ie) ) {


                    for(unsigned int k=kz_b;k<kz_e;k++)
                        for(unsigned int j=jy_b;j<jy_e;j++)
                        {
                            Dxu[IDX_L(ix_e - 2, j, k)] = (-u[IDX_L(ix_e - 3, j, k)]
                                              + u[IDX_L(ix_e - 1, j, k)]
                                             ) * idx_by_2;

                            Dxu[IDX_L(ix_e - 1, j, k)] = (u[IDX_L(ix_e - 3, j, k)]
                                              - 4.0 * u[IDX_L(ix_e - 2, j, k)]
                                              + 3.0 * u[IDX_L(ix_e - 1, j, k)]
                                             ) * idx_by_2;
                        }



                }
    */

#ifdef DEBUG_DERIVS_COMP
            if(isnan(Dxu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
#endif


        }

/*----------------------------------------------------------------------;
 *
 * compute first derivative in y direction
 *
 *----------------------------------------------------------------------*/

    __device__ void _RSWS_deriv42_y(double * const  Dyu, const double * const  u, const double dy, const unsigned int* ijk_lm, const unsigned int * sz, const unsigned int* tile_sz, unsigned int pw, unsigned bflag) {

            const unsigned int i_b=max((int)ijk_lm[2*0+0],(int)3);
            const unsigned int i_e=min((int)ijk_lm[2*0+1],sz[0]-3);

            const unsigned int j_b=ijk_lm[2*1+0]+pw;
            const unsigned int j_e=ijk_lm[2*1+1]-pw;

            const unsigned int k_b=max((int)ijk_lm[2*2+0],int(1));
            const unsigned int k_e=min((int)ijk_lm[2*2+1],sz[2]-1);

            const unsigned int k=(k_b-ijk_lm[4])+threadIdx.y+threadIdx.z*blockDim.y;
            const unsigned int i=(i_b-ijk_lm[0])+threadIdx.x;

            if((i>(i_e-ijk_lm[0])) || (k>(k_e-ijk_lm[4]))) return;

            const double idy = 1.0 / dy;
            const double idy_by_2 = 0.5 * idy;
            const double idy_by_12 = idy / 12.0;


            //printf("dy threadid (%d,%d,%d) loop begin: (%d,%d,%d) loop end: (%d,%d,%d)  tile begin: (%d,%d,%d) tile end: (%d,%d,%d) \n", threadIdx.x,threadIdx.y,threadIdx.z,ix_b,jy_b,kz_b,ix_e,jy_e,kz_e,ijk_lm[0],ijk_lm[2],ijk_lm[4],ijk_lm[1],ijk_lm[3],ijk_lm[5]);

            unsigned int pp = IDX_L(i, (j_b-ijk_lm[2]), k);
            for(unsigned int j=(j_b-ijk_lm[2]);j<=(j_e-ijk_lm[2]);j+=1, pp+=tile_sz[0])
               Dyu[pp] = (u[pp - 2 * tile_sz[0]] - 8.0 * u[pp - tile_sz[0]] + 8.0 * u[pp + tile_sz[0]] - u[pp + 2 * tile_sz[0]]) * idy_by_12;


           /* if ((bflag & (1u << OCT_DIR_DOWN)) && (j_b==jb && (jy_b+ijk_lm[2])==jb) ) {

                    for(unsigned int k=kz_b;k<kz_e;k++)
                        for(unsigned int i=ix_b;i<ix_e;i++)
                        {
                            Dyu[IDX_L(i, jy_b, k)] = (-3.0 * u[IDX_L(i, jy_b, k)]
                                         + 4.0 * u[IDX_L(i, jy_b+1, k)]
                                         - u[IDX_L(i, jy_b+2, k)]
                                        ) * idy_by_2;

                            Dyu[IDX_L(i, jy_b+1, k)] = (-u[IDX_L(i, jy_b, k)]
                                         + u[IDX_L(i, jy_b+2, k)]
                                        ) * idy_by_2;


                        }


                }

            if ((bflag & (1u << OCT_DIR_UP)) && (j_e==je && (jy_e+ijk_lm[2])==je)) {


                    for(unsigned int k=kz_b;k<kz_e;k++)
                        for(unsigned int i=ix_b;i<ix_e;i++)
                        {

                            Dyu[IDX_L(i, jy_e - 2, k)] = (-u[IDX_L(i, jy_e - 3, k)]
                                              + u[IDX_L(i, je - 1, k)]
                                             ) * idy_by_2;

                            Dyu[IDX_L(i, jy_e - 1, k)] = (u[IDX_L(i, jy_e - 3, k)]
                                              - 4.0 * u[IDX_L(i, jy_e - 2, k)]
                                              + 3.0 * u[IDX_L(i, jy_e - 1, k)]
                                             ) * idy_by_2;

                        }



                }
    */

#ifdef DEBUG_DERIVS_COMP
            if(std::isnan(Dyu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
#endif

        }

/*----------------------------------------------------------------------;
 *
 * compute first derivative in z direction
 *
 *----------------------------------------------------------------------*/


    __device__ void _RSWS_deriv42_z(double * const  Dzu, const double * const  u, const double dz, const unsigned int* ijk_lm, const unsigned int * sz, const unsigned int* tile_sz, unsigned int pw, unsigned bflag) {

            const unsigned int i_b=max((int)ijk_lm[2*0+0],(int)3);
            const unsigned int i_e=min((int)ijk_lm[2*0+1],sz[0]-3);

            const unsigned int j_b=max((int)ijk_lm[2*1+0],(int)3);
            const unsigned int j_e=min((int)ijk_lm[2*1+1],sz[1]-3);

            const unsigned int k_b=ijk_lm[2*2+0]+pw;
            const unsigned int k_e=ijk_lm[2*2+1]-pw;

            const unsigned int j=(j_b-ijk_lm[2])+threadIdx.y+threadIdx.z*blockDim.y;
            const unsigned int i=(i_b-ijk_lm[0])+threadIdx.x;

            if((i>(i_e-ijk_lm[0])) || (j>(j_e-ijk_lm[2]))) return;


            const double idz = 1.0 / dz;
            const double idz_by_2 = 0.5 * idz;
            const double idz_by_12 = idz / 12.0;


            //printf("dz threadid (%d,%d,%d) loop begin: (%d,%d,%d) loop end: (%d,%d,%d)  tile begin: (%d,%d,%d) tile end: (%d,%d,%d) \n", threadIdx.x,threadIdx.y,threadIdx.z,ix_b,jy_b,kz_b,ix_e,jy_e,kz_e,ijk_lm[0],ijk_lm[2],ijk_lm[4],ijk_lm[1],ijk_lm[3],ijk_lm[5]);
            unsigned int pp = IDX_L(i, j, (k_b-ijk_lm[4]));
            for(unsigned int k=(k_b-ijk_lm[4]);k<=(k_e-ijk_lm[4]);k+=1, pp+=(tile_sz[0]*tile_sz[1]))
                Dzu[pp] = (u[pp - 2 * tile_sz[0]*tile_sz[1]] - 8.0 * u[pp - tile_sz[0]*tile_sz[1]] + 8.0 * u[pp + tile_sz[0]*tile_sz[1]] - u[pp + 2 * tile_sz[0]*tile_sz[1]]) * idz_by_12;


           /* if ((bflag & (1u << OCT_DIR_BACK)) && (k_b==kb && (kz_b+ijk_lm[4])==kb) ) {

                    for(unsigned int j=jy_b;j<jy_e;j++)
                        for(unsigned int i=ix_b;i<ix_e;i++)
                        {
                            Dzu[IDX_L(i, j, kz_b)] = (-3.0 * u[IDX_L(i, j, kz_b)]
                                         + 4.0 * u[IDX_L(i, j, kz_b+1)]
                                         - u[IDX_L(i, j, 5)]
                                        ) * idz_by_2;

                           Dzu[IDX_L(i, j, kz_b+1)] = (-u[IDX_L(i, j, kz_b)]
                                         + u[IDX_L(i, j, kz_b+2)]
                                        ) * idz_by_2;

                        }

                }

            if ((bflag & (1u << OCT_DIR_FRONT)) && (k_e==ke && (kz_e+ijk_lm[4])==ke) ) {


                    for(unsigned int j=jy_b;j<jy_e;j++)
                        for(unsigned int i=ix_b;i<ix_e;i++)
                        {
                            Dzu[IDX_L(i, j, kz_e - 2)] = (-u[IDX_L(i, j, kz_e - 3)]
                                              + u[IDX_L(i, j, kz_e - 1)]
                                             ) * idz_by_2;

                            Dzu[IDX_L(i, j, kz_e - 1)] = (u[IDX_L(i, j, kz_e - 3)]
                                              - 4.0 * u[IDX_L(i, j, kz_e - 2)]
                                              + 3.0 * u[IDX_L(i, j, kz_e - 1)]
                                             ) * idz_by_2;
                        }



                }*/

#ifdef DEBUG_DERIVS_COMP
            if(std::isnan(Dzu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
#endif


        }


/*----------------------------------------------------------------------;
 *
 * compute second derivative in x direction
 *
 *----------------------------------------------------------------------*/

    __device__ void _RSWS_deriv42_xx(double * const  DxDxu, const double * const  u, const double dx, const unsigned int* ijk_lm, const unsigned int * sz, const unsigned int* tile_sz, unsigned int pw,unsigned bflag) {


            const unsigned int i_b=ijk_lm[2*0+0]+pw;
            const unsigned int i_e=ijk_lm[2*0+1]-pw;

            const unsigned int j_b=max((int)ijk_lm[2*1+0],(int)3);
            const unsigned int j_e=min((int)ijk_lm[2*1+1],sz[1]-3);

            const unsigned int k_b=max((int)ijk_lm[2*2+0],int(3));
            const unsigned int k_e=min((int)ijk_lm[2*2+1],sz[2]-3);

            const unsigned int k=(k_b-ijk_lm[4])+threadIdx.y+threadIdx.z*blockDim.y;
            const unsigned int j=(j_b-ijk_lm[2])+threadIdx.x;

            if((j>(j_e-ijk_lm[2])) || (k>(k_e-ijk_lm[4])) ) return;

            const double idx_sqrd = 1.0 / (dx * dx);
            const double idx_sqrd_by_12 = idx_sqrd / 12.0;

            unsigned int pp = IDX_L((i_b-ijk_lm[0]), j, k);
            for(unsigned int i=(i_b-ijk_lm[0]);i<=(i_e-ijk_lm[0]);i+=1,++pp)
            {


                        DxDxu[pp] = (-u[pp - 2]
                         + 16.0 * u[pp - 1]
                         - 30.0 * u[pp]
                         + 16.0 * u[pp + 1]
                         - u[pp + 2]
                        ) * idx_sqrd_by_12;

            }




           /* if ((bflag & (1u << OCT_DIR_LEFT)) && (i_b==ib && (ix_b+ijk_lm[0])==ib)  ) {

                    for(unsigned int k=kz_b;k<kz_e;k++)
                        for(unsigned int j=jy_b;j<jy_e;j++)
                        {

                            DxDxu[IDX_L(ix_b, j, k)] = (2.0 * u[IDX_L(ix_b, j, k)]
                                           - 5.0 * u[IDX_L(ix_b+1, j, k)]
                                           + 4.0 * u[IDX_L(ix_b+2, j, k)]
                                           - u[IDX_L(ix_b+3, j, k)]
                                          ) * idx_sqrd;

                            DxDxu[IDX_L(ix_b+1, j, k)] = (u[IDX_L(ix_b, j, k)]
                                           - 2.0 * u[IDX_L(ix_b+1, j, k)]
                                           + u[IDX_L(ix_b+2, j, k)]
                                          ) * idx_sqrd;


                        }




                }

            if ((bflag & (1u << OCT_DIR_RIGHT)) && (i_e==ie && (ix_e+ijk_lm[0])==ie) ) {


                    for(unsigned int k=kz_b;k<kz_e;k++)
                        for(unsigned int j=jy_b;j<jy_e;j++)
                        {
                            DxDxu[IDX_L(ix_e - 2, j, k)] = (u[IDX_L(ix_e - 3, j, k)]
                                                - 2.0 * u[IDX_L(ix_e - 2, j, k)]
                                                + u[IDX_L(ix_e - 1, j, k)]
                                               ) * idx_sqrd;

                            DxDxu[IDX_L(ix_e - 1, j, k)] = (-u[IDX_L(ix_e - 4, j, k)]
                                                + 4.0 * u[IDX_L(ix_e - 3, j, k)]
                                                - 5.0 * u[IDX_L(ix_e - 2, j, k)]
                                                + 2.0 * u[IDX_L(ix_e - 1, j, k)]
                                               ) * idx_sqrd;
                        }




                }*/

#ifdef DEBUG_DERIVS_COMP
            if(std::isnan(DxDxu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
#endif

        }


/*----------------------------------------------------------------------;
 *
 * compute second derivative in y direction
 *
 *----------------------------------------------------------------------*/



    __device__ void _RSWS_deriv42_yy(double * const  DyDyu, const double * const  u, const double dy, const unsigned int* ijk_lm, const unsigned int * sz, const unsigned int* tile_sz, unsigned int pw, unsigned bflag) {


            const unsigned int i_b=max((int)ijk_lm[2*0+0],(int)3);
            const unsigned int i_e=min((int)ijk_lm[2*0+1],sz[0]-3);

            const unsigned int j_b=ijk_lm[2*1+0]+pw;
            const unsigned int j_e=ijk_lm[2*1+1]-pw;

            const unsigned int k_b=max((int)ijk_lm[2*2+0],int(3));
            const unsigned int k_e=min((int)ijk_lm[2*2+1],sz[2]-3);

            const unsigned int k=(k_b-ijk_lm[4])+threadIdx.y+threadIdx.z*blockDim.y;
            const unsigned int i=(i_b-ijk_lm[0])+threadIdx.x;

            if((i>(i_e-ijk_lm[0])) || (k>(k_e-ijk_lm[4]))) return;

            const double idy_sqrd = 1.0 / (dy * dy);
            const double idy_sqrd_by_12 = idy_sqrd / 12.0;

            unsigned int pp=IDX_L(i, (j_b-ijk_lm[2]), k);

            for(unsigned int j=(j_b-ijk_lm[2]);j<=(j_e-ijk_lm[2]);j+=1,pp+=tile_sz[0])
            {
                DyDyu[pp] = (-u[pp - 2 * tile_sz[0]] + 16.0 * u[pp - tile_sz[0]] - 30.0 * u[pp]
                            + 16.0 * u[pp + tile_sz[0]] - u[pp + 2 * tile_sz[0]]
                            ) * idy_sqrd_by_12;

            }





           /* if ((bflag & (1u << OCT_DIR_DOWN)) && (j_b==jb && (jy_b+ijk_lm[2])==jb) ) {

                    for(unsigned int k=kz_b;k<kz_e;k++)
                        for(unsigned int i=ix_b;i<ix_e;i++)
                        {
                            DyDyu[IDX_L(i, jy_b, k)] = (2.0 * u[IDX_L(i, jy_b, k)]
                                           - 5.0 * u[IDX_L(i, jy_b+1, k)]
                                           + 4.0 * u[IDX_L(i, jy_b+2, k)]
                                           - u[IDX_L(i, jy_b+3, k)]
                                          ) * idy_sqrd;

                            DyDyu[IDX_L(i, jy_b+1, k)] = (u[IDX_L(i, jy_b, k)]
                                           - 2.0 * u[IDX_L(i, jy_b+1, k)]
                                           + u[IDX_L(i, jy_b+2, k)]
                                          ) * idy_sqrd;

                        }


                }

            if ((bflag & (1u << OCT_DIR_UP)) && (j_e==je && (jy_e+ijk_lm[2])==je)) {

                    for(unsigned int k=kz_b;k<kz_e;k++)
                        for(unsigned int i=ix_b;i<ix_e;i++)
                        {
                            DyDyu[IDX_L(i, jy_e - 2, k)] = (u[IDX_L(i, jy_e - 3, k)]
                                                - 2.0 * u[IDX_L(i, jy_e - 2, k)]
                                                + u[IDX_L(i, jy_e - 1, k)]
                                               ) * idy_sqrd;

                            DyDyu[IDX_L(i, jy_e - 1, k)] = (-u[IDX_L(i, jy_e - 4, k)]
                                                + 4.0 * u[IDX_L(i, jy_e - 3, k)]
                                                - 5.0 * u[IDX_L(i, jy_e - 2, k)]
                                                + 2.0 * u[IDX_L(i, jy_e - 1, k)]
                                               ) * idy_sqrd;

                        }



                }*/

#ifdef DEBUG_DERIVS_COMP
            if(std::isnan(DyDyu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
#endif


        }



/*----------------------------------------------------------------------;
 *
 * compute second derivative in z direction
 *
 *----------------------------------------------------------------------*/


   __device__ void _RSWS_deriv42_zz(double * const  DzDzu, const double * const  u, const double dz, const unsigned int* ijk_lm, const unsigned int * sz, const unsigned int* tile_sz,unsigned int pw, unsigned bflag) {

            const unsigned int i_b=max((int)ijk_lm[2*0+0],(int)3);
            const unsigned int i_e=min((int)ijk_lm[2*0+1],sz[0]-3);

            const unsigned int j_b=max((int)ijk_lm[2*1+0],(int)3);
            const unsigned int j_e=min((int)ijk_lm[2*1+1],sz[1]-3);

            const unsigned int k_b=ijk_lm[2*2+0]+pw;
            const unsigned int k_e=ijk_lm[2*2+1]-pw;

            const unsigned int j=(j_b-ijk_lm[2])+threadIdx.y+threadIdx.z*blockDim.y;
            const unsigned int i=(i_b-ijk_lm[0])+threadIdx.x;

            if((i>(i_e-ijk_lm[0])) || (j>(j_e-ijk_lm[2]))) return;

            const double idz_sqrd = 1.0 / (dz * dz);
            const double idz_sqrd_by_12 = idz_sqrd / 12.0;
            unsigned int pp= IDX_L(i, j, (k_b-ijk_lm[4]));

            for(unsigned int k=(k_b-ijk_lm[4]);k<=(k_e-ijk_lm[4]);k+=1,pp+=(tile_sz[0]*tile_sz[1]))
            {
                DzDzu[pp] = (-u[pp - 2 * tile_sz[0]*tile_sz[1]] + 16.0 * u[pp - tile_sz[0]*tile_sz[1]] - 30.0 * u[pp]
                             + 16.0 * u[pp + tile_sz[0]*tile_sz[1]] - u[pp + 2 * tile_sz[0]*tile_sz[1]]) * idz_sqrd_by_12;
            }



            /*if ((bflag & (1u << OCT_DIR_BACK)) && (k_b==kb && (kz_b+ijk_lm[4])==kb) ) {

                    for(unsigned int j=jy_b;j<jy_e;j++)
                        for(unsigned int i=ix_b;i<ix_e;i++)
                        {
                            DzDzu[IDX_L(i, j, kz_b)] = (2.0 * u[IDX_L(i, j, kz_b)]
                                           - 5.0 * u[IDX_L(i, j, kz_b+1)]
                                           + 4.0 * u[IDX_L(i, j, kz_b+2)]
                                           - u[IDX_L(i, j, kz_b+3)]
                                          ) * idz_sqrd;

                            DzDzu[IDX_L(i, j, kz_b+1)] = (u[IDX_L(i, j, kz_b)]
                                           - 2.0 * u[IDX_L(i, j, kz_b+1)]
                                           + u[IDX_L(i, j, kz_b+2)]
                                          ) * idz_sqrd;

                        }


                }

            if ((bflag & (1u << OCT_DIR_FRONT)) && (k_e==ke && (kz_e+ijk_lm[4])==ke) ) {


                    for(unsigned int j=jy_b;j<jy_e;j++)
                        for(unsigned int i=ix_b;i<ix_e;i++)
                        {
                            DzDzu[IDX_L(i, j, kz_e - 2)] = (u[IDX_L(i, j, kz_e - 3)]
                                                - 2.0 * u[IDX_L(i, j, kz_e - 2)]
                                                + u[IDX_L(i, j, kz_e - 1)]
                                               ) * idz_sqrd;

                            DzDzu[IDX_L(i, j, kz_e - 1)] = (-u[IDX_L(i, j, kz_e - 4)]
                                                + 4.0 * u[IDX_L(i, j, kz_e - 3)]
                                                - 5.0 * u[IDX_L(i, j, kz_e - 2)]
                                                + 2.0 * u[IDX_L(i, j, kz_e - 1)]
                                               ) * idz_sqrd;

                        }


                }*/

#ifdef DEBUG_DERIVS_COMP
            if(std::isnan(DzDzu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
#endif


        }


/*----------------------------------------------------------------------;
 *
 * compute first advective derivative in x direction
 *
 *----------------------------------------------------------------------*/
    __device__    void _RSWS_deriv42adv_x(double * const  Dxu, const double * const  u, const double dx, const unsigned int* ijk_lm, const unsigned int * sz, const unsigned int* tile_sz, const bool * const betax, unsigned int pw, unsigned bflag) {

            const unsigned int i_b=ijk_lm[2*0+0]+pw;
            const unsigned int i_e=ijk_lm[2*0+1]-pw;

            const unsigned int j_b=max((int)ijk_lm[2*1+0],(int)3);
            const unsigned int j_e=min((int)ijk_lm[2*1+1],sz[1]-3);

            const unsigned int k_b=max((int)ijk_lm[2*2+0],int(3));
            const unsigned int k_e=min((int)ijk_lm[2*2+1],sz[2]-3);

            const unsigned int k=(k_b-ijk_lm[4])+threadIdx.y+threadIdx.z*blockDim.y;
            const unsigned int j=(j_b-ijk_lm[2])+threadIdx.x;

            if((j>(j_e-ijk_lm[2])) || (k>(k_e-ijk_lm[4])) ) return;

            const double idx = 1.0 / dx;
            const double idx_by_2 = 0.50 * idx;
            const double idx_by_12 = idx / 12.0;
            unsigned int pp=IDX_L((i_b-ijk_lm[0]), j, k);

            for(unsigned int i=(i_b-ijk_lm[0]);i<=(i_e-ijk_lm[0]);i+=1,++pp)
            {
                        (betax[pp]) ? Dxu[pp] = (-3.0 * u[pp - 1]
                                                 - 10.0 * u[pp]
                                                 + 18.0 * u[pp + 1]
                                                 - 6.0 * u[pp + 2]
                                                 + u[pp + 3]
                                                ) * idx_by_12 : Dxu[pp] = (-u[pp - 3]
                                                                           + 6.0 * u[pp - 2]
                                                                           - 18.0 * u[pp - 1]
                                                                           + 10.0 * u[pp]
                                                                           + 3.0 * u[pp + 1]
                                                                           ) * idx_by_12;

            }



            /*if ((bflag & (1u << OCT_DIR_LEFT)) && (i_b==ib && (ix_b+ijk_lm[0])==ib)  ) {


                    for(unsigned int k=kz_b;k<kz_e;k++)
                        for(unsigned int j=jy_b;j<jy_e;j++)
                        {
                            Dxu[IDX_L(ix_b, j, k)] = (-3.0 * u[IDX_L(ix_b, j, k)]
                                         + 4.0 * u[IDX_L(ix_b+1, j, k)]
                                         - u[IDX_L(ix_b+2, j, k)]
                                        ) * idx_by_2;

                            if (betax[IDX_L(ix_b+1, j, k)] > 0.0) {
                                Dxu[IDX_L(ix_b+1, j, k)] = (-3.0 * u[IDX_L(ix_b+1, j, k)]
                                                    + 4.0 * u[IDX_L(ix_b+2, j, k)]
                                                    - u[IDX_L(ix_b+3, j, k)]
                                                    ) * idx_by_2;
                            } else {
                                Dxu[IDX_L(ix_b+1, j, k)] = (-u[IDX_L(ix_b, j, k)]
                                                    + u[IDX_L(ix_b+2, j, k)]
                                                    ) * idx_by_2;
                            }

                            if (betax[IDX_L(ix_b+2, j, k)] > 0.0) {
                                Dxu[IDX_L(ix_b+2, j, k)] = (-3.0 * u[IDX_L(ix_b+1, j, k)]
                                                    - 10.0 * u[IDX_L(ix_b+2, j, k)]
                                                    + 18.0 * u[IDX_L(ix_b+3, j, k)]
                                                    - 6.0 * u[IDX_L(ix_b+4, j, k)]
                                                    + u[IDX_L(ix_b+5, j, k)]
                                                    ) * idx_by_12;
                            } else {
                                Dxu[IDX_L(ix_b+2, j, k)] = (u[IDX_L(ix_b, j, k)]
                                                    - 4.0 * u[IDX_L(ix_b+1, j, k)]
                                                    + 3.0 * u[IDX_L(ix_b+2, j, k)]
                                                    ) * idx_by_2;
                            }

                        }




                }

            if ((bflag & (1u << OCT_DIR_RIGHT)) && (i_e==ie && (ix_e+ijk_lm[0])==ie) ) {


                    for(unsigned int k=kz_b;k<kz_e;k++)
                        for(unsigned int j=jy_b;j<jy_e;j++)
                        {
                            if (betax[IDX_L(ix_e - 3, j, k)] < 0.0) {
                                Dxu[IDX_L(ix_e - 3, j, k)] = (-3.0 * u[IDX_L(ix_e - 3, j, k)]
                                                          + 4.0 * u[IDX_L(ix_e - 2, j, k)]
                                                          - u[IDX_L(ix_e - 1, j, k)]
                                                         ) * idx_by_2;
                            } else {
                                Dxu[IDX_L(ix_e - 3, j, k)] = (-u[IDX_L(ix_e - 6, j, k)]
                                                          + 6.0 * u[IDX_L(ix_e - 5, j, k)]
                                                          - 18.0 * u[IDX_L(ix_e - 4, j, k)]
                                                          + 10.0 * u[IDX_L(ix_e - 3, j, k)]
                                                          + 3.0 * u[IDX_L(ix_e - 2, j, k)]
                                                         ) * idx_by_12;
                            }

                            if (betax[IDX_L(ix_e - 2, j, k)] > 0.0) {
                                Dxu[IDX_L(ix_e - 2, j, k)] = (-u[IDX_L(ix_e - 3, j, k)]
                                                          + u[IDX_L(ix_e - 1, j, k)]
                                                         ) * idx_by_2;
                            } else {
                                Dxu[IDX_L(ix_e - 2, j, k)] = (u[IDX_L(ix_e - 4, j, k)]
                                                          - 4.0 * u[IDX_L(ix_e - 3, j, k)]
                                                          + 3.0 * u[IDX_L(ix_e - 2, j, k)]
                                                         ) * idx_by_2;
                            }

                            Dxu[IDX_L(ix_e - 1, j, k)] = (u[IDX_L(ix_e - 3, j, k)]
                                                      - 4.0 * u[IDX_L(ix_e - 2, j, k)]
                                                      + 3.0 * u[IDX_L(ix_e - 1, j, k)]
                                                     ) * idx_by_2;

                        }




                }*/

#ifdef DEBUG_DERIVS_COMP
            if(std::isnan(Dxu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
#endif


        }


/*----------------------------------------------------------------------;
 *
 * compute first advective derivative in y direction
 *
 *----------------------------------------------------------------------*/
    __device__  void _RSWS_deriv42adv_y(double * const  Dyu, const double * const  u, const double dy, const unsigned int* ijk_lm, const unsigned int * sz, const unsigned int* tile_sz, const bool * const betay, unsigned int pw, unsigned bflag) {

            const unsigned int i_b=max((int)ijk_lm[2*0+0],(int)3);
            const unsigned int i_e=min((int)ijk_lm[2*0+1],sz[0]-3);

            const unsigned int j_b=ijk_lm[2*1+0]+pw;
            const unsigned int j_e=ijk_lm[2*1+1]-pw;

            const unsigned int k_b=max((int)ijk_lm[2*2+0],int(3));
            const unsigned int k_e=min((int)ijk_lm[2*2+1],sz[2]-3);

            const unsigned int k=(k_b-ijk_lm[4])+threadIdx.y+threadIdx.z*blockDim.y;
            const unsigned int i=(i_b-ijk_lm[0])+threadIdx.x;

            if((i>(i_e-ijk_lm[0])) || (k>(k_e-ijk_lm[4]))) return;

            const double idy = 1.0 / dy;
            const double idy_by_2 = 0.50 * idy;
            const double idy_by_12 = idy / 12.0;

            unsigned int pp=IDX_L(i, (j_b-ijk_lm[2]), k);

            for(unsigned int j=(j_b-ijk_lm[2]);j<=(j_e-ijk_lm[2]);j+=1,pp+=(tile_sz[0]))
            {
                         (betay[pp]) ? Dyu[pp] = (-3.0 * u[pp - tile_sz[0]]
                                    - 10.0 * u[pp]
                                    + 18.0 * u[pp + tile_sz[0]]
                                    - 6.0 * u[pp + 2 * tile_sz[0]]
                                    + u[pp + 3 * tile_sz[0]]
                                    ) * idy_by_12:
                                            Dyu[pp] = (-u[pp - 3 * tile_sz[0]]
                                            + 6.0 * u[pp - 2 * tile_sz[0]]
                                            - 18.0 * u[pp - tile_sz[0]]
                                            + 10.0 * u[pp]
                                            + 3.0 * u[pp + tile_sz[0]]
                                            ) * idy_by_12;

            }


            /*if ((bflag & (1u << OCT_DIR_DOWN)) && (j_b==jb && (jy_b+ijk_lm[2])==jb) ) {

                    for(unsigned int k=kz_b;k<kz_e;k++)
                        for(unsigned int i=ix_b;i<ix_e;i++)
                        {
                            Dyu[IDX_L(i, jy_b, k)] = (-3.0 * u[IDX_L(i, jy_b, k)]
                                         + 4.0 * u[IDX_L(i, jy_b+1, k)]
                                         - u[IDX_L(i, jy_b+2, k)]
                                        ) * idy_by_2;

                            if (betay[IDX_L(i, jy_b+1, k)] > 0.0) {
                                Dyu[IDX_L(i, jy_b+1, k)] = (-3.0 * u[IDX_L(i, jy_b+1, k)]
                                                    + 4.0 * u[IDX_L(i, jy_b+2, k)]
                                                    - u[IDX_L(i, jy_b+3, k)]
                                                    ) * idy_by_2;
                            } else {
                                Dyu[IDX_L(i, jy_b+1, k)] = (-u[IDX_L(i, jy_b, k)]
                                                    + u[IDX_L(i, jy_b+2, k)]
                                                    ) * idy_by_2;
                            }

                            if (betay[IDX_L(i, jy_b+2, k)] > 0.0) {
                                Dyu[IDX_L(i, jy_b+2, k)] = (-3.0 * u[IDX_L(i, jy_b+1, k)]
                                                    - 10.0 * u[IDX_L(i, jy_b+2, k)]
                                                    + 18.0 * u[IDX_L(i, jy_b+3, k)]
                                                    - 6.0 * u[IDX_L(i, jy_b+4, k)]
                                                    + u[IDX_L(i, jy_b+5, k)]
                                                    ) * idy_by_12;
                            } else {
                                Dyu[IDX_L(i, jy_b+2, k)] = (u[IDX_L(i, jy_b, k)]
                                                    - 4.0 * u[IDX_L(i, jy_b+1, k)]
                                                    + 3.0 * u[IDX_L(i, jy_b+2, k)]
                                                    ) * idy_by_2;
                            }

                        }



                }

            if ((bflag & (1u << OCT_DIR_UP)) && (j_e==je && (jy_e+ijk_lm[2])==je)) {


                    for(unsigned int k=kz_b;k<kz_e;k++)
                        for(unsigned int i=ix_b;i<ix_e;i++)
                        {
                            if (betay[IDX_L(i, jy_e - 3, k)] < 0.0) {
                                Dyu[IDX_L(i, jy_e - 3, k)] = (-3.0 * u[IDX_L(i, jy_e - 3, k)]
                                                          + 4.0 * u[IDX_L(i, jy_e - 2, k)]
                                                          - u[IDX_L(i, jy_e - 1, k)]
                                                         ) * idy_by_2;
                            } else {
                                Dyu[IDX_L(i, jy_e - 3, k)] = (-u[IDX_L(i, jy_e - 6, k)]
                                                          + 6.0 * u[IDX_L(i, jy_e - 5, k)]
                                                          - 18.0 * u[IDX_L(i, jy_e - 4, k)]
                                                          + 10.0 * u[IDX_L(i, jy_e - 3, k)]
                                                          + 3.0 * u[IDX_L(i, jy_e - 2, k)]
                                                         ) * idy_by_12;
                            }

                            if (betay[IDX_L(i, jy_e - 2, k)] > 0.0) {
                                Dyu[IDX_L(i, jy_e - 2, k)] = (-u[IDX_L(i, jy_e - 3, k)]
                                                          + u[IDX_L(i, jy_e - 1, k)]
                                                         ) * idy_by_2;
                            } else {
                                Dyu[IDX_L(i, jy_e - 2, k)] = (u[IDX_L(i, jy_e - 4, k)]
                                                          - 4.0 * u[IDX_L(i, jy_e - 3, k)]
                                                          + 3.0 * u[IDX_L(i, jy_e - 2, k)]
                                                         ) * idy_by_2;
                            }

                            Dyu[IDX_L(i, jy_e - 1, k)] = (u[IDX_L(i, jy_e - 3, k)]
                                                      - 4.0 * u[IDX_L(i, jy_e - 2, k)]
                                                      + 3.0 * u[IDX_L(i, jy_e - 1, k)]
                                                     ) * idy_by_2;

                        }




                }*/


#ifdef DEBUG_DERIVS_COMP
            if(std::isnan(Dyu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
#endif


        }


/*----------------------------------------------------------------------;
 *
 * compute first advective derivative in z direction
 *
 *----------------------------------------------------------------------*/


    __device__  void _RSWS_deriv42adv_z(double * const  Dzu, const double * const  u, const double dz, const unsigned int* ijk_lm, const unsigned int * sz, const unsigned int* tile_sz, const bool * const betaz, unsigned int pw, unsigned bflag) {


            const unsigned int i_b=max((int)ijk_lm[2*0+0],(int)3);
            const unsigned int i_e=min((int)ijk_lm[2*0+1],sz[0]-3);

            const unsigned int j_b=max((int)ijk_lm[2*1+0],(int)3);
            const unsigned int j_e=min((int)ijk_lm[2*1+1],sz[1]-3);

            const unsigned int k_b=ijk_lm[2*2+0]+pw;
            const unsigned int k_e=ijk_lm[2*2+1]-pw;

            const unsigned int j=(j_b-ijk_lm[2])+threadIdx.y+threadIdx.z*blockDim.y;
            const unsigned int i=(i_b-ijk_lm[0])+threadIdx.x;

            if((i>(i_e-ijk_lm[0])) || (j>(j_e-ijk_lm[2]))) return;

            const double idz = 1.0 / dz;
            const double idz_by_2 = 0.50 * idz;
            const double idz_by_12 = idz / 12.0;

            unsigned int pp = IDX_L(i, j, (k_b-ijk_lm[4]));

            for(unsigned int k=(k_b-ijk_lm[4]);k<=(k_e-ijk_lm[4]);k+=1,pp+=(tile_sz[0]*tile_sz[1]))
            {
                    (betaz[pp])? Dzu[pp]=(-3.0 * u[pp - tile_sz[0]*tile_sz[1]]
                                       - 10.0 * u[pp]
                                       + 18.0 * u[pp + tile_sz[0]*tile_sz[1]]
                                       - 6.0 * u[pp + 2 * tile_sz[0]*tile_sz[1]]
                                       + u[pp + 3 * tile_sz[0]*tile_sz[1]]
                                      ) * idz_by_12 : Dzu[pp] = (-u[pp - 3 * tile_sz[0]*tile_sz[1]]
                                                                 + 6.0 * u[pp - 2 * tile_sz[0]*tile_sz[1]]
                                                                 - 18.0 * u[pp - tile_sz[0]*tile_sz[1]]
                                                                 + 10.0 * u[pp]
                                                                 + 3.0 * u[pp + tile_sz[0]*tile_sz[1]] ) * idz_by_12;


            }





            /*if ((bflag & (1u << OCT_DIR_BACK)) && (k_b==kb && (kz_b+ijk_lm[4])==kb) ) {

                    for(unsigned int j=jy_b;j<jy_e;j++)
                        for(unsigned int i=ix_b;i<ix_e;i++)
                        {

                            Dzu[IDX_L(i, j, kz_b)] = (-3.0 * u[IDX_L(i, j, kz_b)]
                                         + 4.0 * u[IDX_L(i, j, kz_b+1)]
                                         - u[IDX_L(i, j, kz_b+2)]
                                        ) * idz_by_2;

                            if (betaz[IDX_L(i, j, kz_b+1)] > 0.0) {
                                Dzu[IDX_L(i, j, kz_b+1)] = (-3.0 * u[IDX_L(i, j, kz_b+1)]
                                                    + 4.0 * u[IDX_L(i, j, kz_b+2)]
                                                    - u[IDX_L(i, j, kz_b+3)]
                                                    ) * idz_by_2;
                            } else {
                                Dzu[IDX_L(i, j, kz_b+1)] = (-u[IDX_L(i, j, kz_b)]
                                                    + u[IDX_L(i, j, kz_b+2)]
                                                    ) * idz_by_2;
                            }

                            if (betaz[IDX_L(i, j, kz_b+2)] > 0.0) {
                                Dzu[IDX_L(i, j, kz_b+2)] = (-3.0 * u[IDX_L(i, j, kz_b+1)]
                                                    - 10.0 * u[IDX_L(i, j, kz_b+2)]
                                                    + 18.0 * u[IDX_L(i, j, kz_b+3)]
                                                    - 6.0 * u[IDX_L(i, j, kz_b+4)]
                                                    + u[IDX_L(i, j, kz_b+5)]
                                                    ) * idz_by_12;
                            } else {
                                Dzu[IDX_L(i, j, kz_b+2)] = (u[IDX_L(i, j, kz_b)]
                                                    - 4.0 * u[IDX_L(i, j, kz_b+1)]
                                                    + 3.0 * u[IDX_L(i, j, kz_b+2)]
                                                    ) * idz_by_2;
                            }

                        }



                }

            if ((bflag & (1u << OCT_DIR_FRONT)) && (k_e==ke && (kz_e+ijk_lm[4])==ke) ) {

                    for(unsigned int j=jy_b;j<jy_e;j++)
                        for(unsigned int i=ix_b;i<ix_e;i++)
                        {

                            if (betaz[IDX_L(i, j, kz_e - 3)] < 0.0) {
                                Dzu[IDX_L(i, j, kz_e - 3)] = (-3.0 * u[IDX_L(i, j, kz_e - 3)]
                                                          + 4.0 * u[IDX_L(i, j, kz_e - 2)]
                                                          - u[IDX_L(i, j, kz_e - 1)]
                                                         ) * idz_by_2;
                            } else {
                                Dzu[IDX_L(i, j, kz_e - 3)] = (-u[IDX_L(i, j, kz_e - 6)]
                                                          + 6.0 * u[IDX_L(i, j, kz_e - 5)]
                                                          - 18.0 * u[IDX_L(i, j, kz_e - 4)]
                                                          + 10.0 * u[IDX_L(i, j, kz_e - 3)]
                                                          + 3.0 * u[IDX_L(i, j, kz_e - 2)]
                                                         ) * idz_by_12;
                            }

                            if (betaz[IDX_L(i, j, kz_e - 2)] > 0.0) {
                                Dzu[IDX_L(i, j, kz_e - 2)] = (-u[IDX_L(i, j, kz_e - 3)]
                                                          + u[IDX_L(i, j, kz_e - 1)]
                                                         ) * idz_by_2;
                            } else {
                                Dzu[IDX_L(i, j, kz_e - 2)] = (u[IDX_L(i, j, kz_e - 4)]
                                                          - 4.0 * u[IDX_L(i, j, kz_e - 3)]
                                                          + 3.0 * u[IDX_L(i, j, kz_e - 2)]
                                                         ) * idz_by_2;
                            }

                            Dzu[IDX_L(i, j, ke - 1)] = (u[IDX_L(i, j, ke - 3)]
                                                      - 4.0 * u[IDX_L(i, j, ke - 2)]
                                                      + 3.0 * u[IDX_L(i, j, ke - 1)]
                                                     ) * idz_by_2;

                        }

                }*/

#ifdef DEBUG_DERIVS_COMP
            if(std::isnan(Dzu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
#endif


        }

/*----------------------------------------------------------------------
 *
 * compute Kriess-Oliger derivative in x direction
 *
 *----------------------------------------------------------------------*/


    __device__  void _RSWS_ko_deriv42_x(double * const  Du, const double * const  u, const double dx, const unsigned int* ijk_lm, const unsigned int * sz, const unsigned int* tile_sz, unsigned int pw, unsigned bflag) {


            const unsigned int i_b=ijk_lm[2*0+0]+pw;
            const unsigned int i_e=ijk_lm[2*0+1]-pw;

            const unsigned int j_b=max((int)ijk_lm[2*1+0],(int)3);
            const unsigned int j_e=min((int)ijk_lm[2*1+1],sz[1]-3);

            const unsigned int k_b=max((int)ijk_lm[2*2+0],int(3));
            const unsigned int k_e=min((int)ijk_lm[2*2+1],sz[2]-3);

            const unsigned int k=(k_b-ijk_lm[4])+threadIdx.y+threadIdx.z*blockDim.y;
            const unsigned int j=(j_b-ijk_lm[2])+threadIdx.x;

            if((j>(j_e-ijk_lm[2])) || (k>(k_e-ijk_lm[4])) ) return;


            const double pre_factor_6_dx = -1.0 / 64.0 / dx;

            const double smr3 = 59.0 / 48.0 * 64 * dx;
            const double smr2 = 43.0 / 48.0 * 64 * dx;
            const double smr1 = 49.0 / 48.0 * 64 * dx;

            unsigned int pp = IDX_L((i_b-ijk_lm[0]), j, k);

            /*const double spr3 = smr3;
            const double spr2 = smr2;
            const double spr1 = smr1;*/

            for(unsigned int i=(i_b-ijk_lm[0]);i<=(i_e-ijk_lm[0]);i+=1,++pp)
            {
                  Du[pp] = pre_factor_6_dx *
                       (
                                         -u[pp - 3]
                                         + 6.0 * u[pp - 2]
                                         - 15.0 * u[pp - 1]
                                         + 20.0 * u[pp]
                                         - 15.0 * u[pp + 1]
                                         + 6.0 * u[pp + 2]
                                         - u[pp + 3]
                       );
            }



            if(ijk_lm[0]==pw)
            {


                 Du[IDX_L((i_b-ijk_lm[0]), j, k)] = pre_factor_6_dx *
                                          (
                                                -u[IDX_L((i_b-ijk_lm[0]) + 4, j, k)]
                                                + 6.0 * u[IDX_L((i_b-ijk_lm[0]) + 3, j, k)]
                                                - 15.0 * u[IDX_L((i_b-ijk_lm[0]) + 2, j, k)]
                                                + 20.0 * u[IDX_L((i_b-ijk_lm[0]) + 1, j, k)]
                                                - 15.0 * u[IDX_L((i_b-ijk_lm[0]), j, k)]
                                                + 6.0 * u[IDX_L((i_b-ijk_lm[0]) - 1, j, k)]
                                                - u[IDX_L((i_b-ijk_lm[0]) - 2, j, k)]
                                          );


            }


            if(ijk_lm[1]==(sz[0]-1-pw))
            {

                 Du[IDX_L((i_e-ijk_lm[0]), j, k)] = pre_factor_6_dx *
                                             (
                                                     -u[IDX_L((i_e-ijk_lm[0]) + 1, j, k)]
                                                     + 6.0 * u[IDX_L((i_e-ijk_lm[0]), j, k)]
                                                     - 15.0 * u[IDX_L((i_e-ijk_lm[0]) - 1, j, k)]
                                                     + 20.0 * u[IDX_L((i_e-ijk_lm[0]) - 2, j, k)]
                                                     - 15.0 * u[IDX_L((i_e-ijk_lm[0]) - 3, j, k)]
                                                     + 6.0 * u[IDX_L((i_e-ijk_lm[0]) - 4, j, k)]
                                                     - u[IDX_L((i_e-ijk_lm[0]) - 5, j, k)]
                                             );

            }


        /*if ((bflag & (1u << OCT_DIR_LEFT)) && (i_b==ib && (ix_b+ijk_lm[0])==ib)  ) {


                for(unsigned int k=kz_b;k<kz_e;k++)
                    for(unsigned int j=jy_b;j<jy_e;j++)
                    {
                        Du[IDX_L(ix_b, j, k)] = (u[IDX_L(ix_b+3, j, k)]
                                    - 3.0 * u[IDX_L(ix_b+2, j, k)]
                                    + 3.0 * u[IDX_L(ix_b+1, j, k)]
                                    - u[IDX_L(ix_b, j, k)]
                                   ) / smr3;

                        Du[IDX_L(ix_b+1, j, k)] = (
                                                u[IDX_L(ix_b+4, j, k)]
                                                - 6.0 * u[IDX_L(ix_b+3, j, k)]
                                                + 12.0 * u[IDX_L(ix_b+2, j, k)]
                                                - 10.0 * u[IDX_L(ix_b+1, j, k)]
                                                + 3.0 * u[IDX_L(ix_b, j, k)]
                                        ) / smr2;
                        Du[IDX_L(ix_b+2, j, k)] = (
                                                u[IDX_L(ix_b+5, j, k)]
                                                - 6.0 * u[IDX_L(ix_b+4, j, k)]
                                                + 15.0 * u[IDX_L(ix_b+3, j, k)]
                                                - 19.0 * u[IDX_L(ix_b+2, j, k)]
                                                + 12.0 * u[IDX_L(ix_b+1, j, k)]
                                                - 3.0 * u[IDX_L(ix_b, j, k)]
                                        ) / smr1;
                            }




            }

        if ((bflag & (1u << OCT_DIR_RIGHT)) && (i_e==ie && (ix_e+ijk_lm[0])==ie) ) {

                for(unsigned int k=kz_b;k<kz_e;k++)
                    for(unsigned int j=jy_b;j<jy_e;j++)
                    {
                        Du[IDX_L(ix_e - 3, j, k)] = (
                            u[IDX_L(ix_e - 6, j, k)]
                            - 6.0 * u[IDX_L(ix_e - 5, j, k)]
                            + 15.0 * u[IDX_L(ix_e - 4, j, k)]
                            - 19.0 * u[IDX_L(ix_e - 3, j, k)]
                            + 12.0 * u[IDX_L(ix_e - 2, j, k)]
                            - 3.0 * u[IDX_L(ix_e - 1, j, k)]
                    ) / spr1;

                    Du[IDX_L(ix_e - 2, j, k)] = (
                                                u[IDX_L(ix_e - 5, j, k)]
                                                - 6.0 * u[IDX_L(ix_e - 4, j, k)]
                                                + 12.0 * u[IDX_L(ix_e - 3, j, k)]
                                                - 10.0 * u[IDX_L(ix_e - 2, j, k)]
                                                + 3.0 * u[IDX_L(ix_e - 1, j, k)]
                                        ) / spr2;

                    Du[IDX_L(ix_e - 1, j, k)] = (
                                                u[IDX_L(ix_e - 4, j, k)]
                                                - 3.0 * u[IDX_L(ix_e - 3, j, k)]
                                                + 3.0 * u[IDX_L(ix_e - 2, j, k)]
                                                - u[IDX_L(ix_e - 1, j, k)]
                                        ) / spr3;
                    }



            }*/

#ifdef DEBUG_DERIVS_COMP
            if(std::isnan(Du[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
#endif


        }



/*----------------------------------------------------------------------
 *
 * compute Kriess-Oliger derivative in y direction
 *
 *----------------------------------------------------------------------*/

    __device__  void _RSWS_ko_deriv42_y(double * const  Du, const double * const  u, const double dy, const unsigned int* ijk_lm, const unsigned int * sz, const unsigned int* tile_sz, unsigned int pw, unsigned bflag) {



            const unsigned int i_b=max((int)ijk_lm[2*0+0],(int)3);
            const unsigned int i_e=min((int)ijk_lm[2*0+1],sz[0]-3);

            const unsigned int j_b=ijk_lm[2*1+0]+pw;
            const unsigned int j_e=ijk_lm[2*1+1]-pw;

            const unsigned int k_b=max((int)ijk_lm[2*2+0],int(3));
            const unsigned int k_e=min((int)ijk_lm[2*2+1],sz[2]-3);

            const unsigned int k=(k_b-ijk_lm[4])+threadIdx.y+threadIdx.z*blockDim.y;
            const unsigned int i=(i_b-ijk_lm[0])+threadIdx.x;

            if((i>(i_e-ijk_lm[0])) || (k>(k_e-ijk_lm[4]))) return;

            const double pre_factor_6_dy = -1.0 / 64.0 / dy;

            const double smr3 = 59.0 / 48.0 * 64 * dy;
            const double smr2 = 43.0 / 48.0 * 64 * dy;
            const double smr1 = 49.0 / 48.0 * 64 * dy;

            /*const double spr3 = smr3;
            const double spr2 = smr2;
            const double spr1 = smr1;*/

            unsigned int pp=IDX_L(i, (j_b-ijk_lm[2]), k);
            for(unsigned int j=(j_b-ijk_lm[2]);j<=(j_e-ijk_lm[2]);j+=1,pp+=(tile_sz[0]*tile_sz[1]))
            {
                        Du[pp] = pre_factor_6_dy *
                                (
                                        -u[pp - 3 * tile_sz[0]]
                                        + 6.0 * u[pp - 2 * tile_sz[0]]
                                        - 15.0 * u[pp - tile_sz[0]]
                                        + 20.0 * u[pp]
                                        - 15.0 * u[pp + tile_sz[0]]
                                        + 6.0 * u[pp + 2 * tile_sz[0]]
                                        - u[pp + 3 * tile_sz[0]]
                                );

            }


            if((ijk_lm[2])==pw)
            {


                        Du[IDX_L(i, (j_b-ijk_lm[2]), k)] = pre_factor_6_dy *
                                              (
                                                      -u[IDX_L(i, (j_b-ijk_lm[2]) + 4, k)]
                                                      + 6.0 * u[IDX_L(i, (j_b-ijk_lm[2]) + 3, k)]
                                                      - 15.0 * u[IDX_L(i, (j_b-ijk_lm[2]) + 2, k)]
                                                      + 20.0 * u[IDX_L(i, (j_b-ijk_lm[2]) + 1, k)]
                                                      - 15.0 * u[IDX_L(i, (j_b-ijk_lm[2]), k)]
                                                      + 6.0 * u[IDX_L(i, (j_b-ijk_lm[2]) - 1, k)]
                                                      - u[IDX_L(i, (j_b-ijk_lm[2]) - 2, k)]
                                              );

            }


            if(ijk_lm[3]==(sz[1]-1-pw))
            {

                        Du[IDX_L(i, (j_e-ijk_lm[2]), k)] = pre_factor_6_dy *
                                                  (
                                                          -u[IDX_L(i, (j_e-ijk_lm[2]) + 1, k)]
                                                          + 6.0 * u[IDX_L(i, (j_e-ijk_lm[2]), k)]
                                                          - 15.0 * u[IDX_L(i, (j_e-ijk_lm[2]) - 1, k)]
                                                          + 20.0 * u[IDX_L(i, (j_e-ijk_lm[2]) - 2, k)]
                                                          - 15.0 * u[IDX_L(i, (j_e-ijk_lm[2]) - 3, k)]
                                                          + 6.0 * u[IDX_L(i, (j_e-ijk_lm[2]) - 4, k)]
                                                          - u[IDX_L(i, (j_e-ijk_lm[2]) - 5, k)]
                                                  );

            }


        /*if ((bflag & (1u << OCT_DIR_DOWN)) && (j_b==jb && (jy_b+ijk_lm[2])==jb) ) {

                for(unsigned int k=kz_b;k<kz_e;k++)
                    for(unsigned int i=ix_b;i<ix_e;i++)
                    {
                        Du[IDX_L(i, jy_b, k)] = (u[IDX_L(i, jy_b+3, k)]
                        - 3.0 * u[IDX_L(i, jy_b+2, k)]
                        + 3.0 * u[IDX_L(i, jy_b+1, k)]
                        - u[IDX_L(i, jy_b, k)]
                       ) / smr3;

                       Du[IDX_L(i, jy_b+1, k)] = (
                                                u[IDX_L(i, jy_b+4, k)]
                                                - 6.0 * u[IDX_L(i, jy_b+3, k)]
                                                + 12.0 * u[IDX_L(i, jy_b+2, k)]
                                                - 10.0 * u[IDX_L(i, jy_b+1, k)]
                                                + 3.0 * u[IDX_L(i, jy_b, k)]
                                        ) / smr2;
                       Du[IDX_L(i, jy_b+2, k)] = (
                                                u[IDX_L(i, jy_b+5, k)]
                                                - 6.0 * u[IDX_L(i, jy_b+4, k)]
                                                + 15.0 * u[IDX_L(i, jy_b+3, k)]
                                                - 19.0 * u[IDX_L(i, jy_b+2, k)]
                                                + 12.0 * u[IDX_L(i, jy_b+1, k)]
                                                - 3.0 * u[IDX_L(i, jy_b, k)]
                                        ) / smr1;

                    }



            }

        if ((bflag & (1u << OCT_DIR_UP)) && (j_e==je && (jy_e+ijk_lm[2])==je)) {

                for(unsigned int k=kz_b;k<kz_e;k++)
                    for(unsigned int i=ix_b;i<ix_e;i++)
                    {
                        Du[IDX_L(i, jy_e - 3, k)] = (
                            u[IDX_L(i, jy_e - 6, k)]
                            - 6.0 * u[IDX_L(i, jy_e - 5, k)]
                            + 15.0 * u[IDX_L(i, jy_e - 4, k)]
                            - 19.0 * u[IDX_L(i, jy_e - 3, k)]
                            + 12.0 * u[IDX_L(i, jy_e - 2, k)]
                            - 3.0 * u[IDX_L(i, jy_e - 1, k)]
                    ) / spr1;

                    Du[IDX_L(i, jy_e - 2, k)] = (
                                                u[IDX_L(i, jy_e - 5, k)]
                                                - 6.0 * u[IDX_L(i, jy_e - 4, k)]
                                                + 12.0 * u[IDX_L(i, jy_e - 3, k)]
                                                - 10.0 * u[IDX_L(i, jy_e - 2, k)]
                                                + 3.0 * u[IDX_L(i, jy_e - 1, k)]
                                        ) / spr2;

                    Du[IDX_L(i, jy_e - 1, k)] = (
                                                u[IDX_L(i, jy_e - 4, k)]
                                                - 3.0 * u[IDX_L(i, jy_e - 3, k)]
                                                + 3.0 * u[IDX_L(i, jy_e - 2, k)]
                                                - u[IDX_L(i, jy_e - 1, k)]
                                        ) / spr3;

                    }



            }*/

#ifdef DEBUG_DERIVS_COMP
            if(std::isnan(Du[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
#endif


        }



/*----------------------------------------------------------------------
 *
 * compute Kriess-Oliger derivative in z direction
 *
 *----------------------------------------------------------------------*/

    __device__  void _RSWS_ko_deriv42_z(double * const  Du, const double * const  u, const double dz, const unsigned int* ijk_lm, const unsigned int * sz, const unsigned int* tile_sz, unsigned int pw, unsigned bflag) {



            const unsigned int i_b=max((int)ijk_lm[2*0+0],(int)3);
            const unsigned int i_e=min((int)ijk_lm[2*0+1],sz[0]-3);

            const unsigned int j_b=max((int)ijk_lm[2*1+0],(int)3);
            const unsigned int j_e=min((int)ijk_lm[2*1+1],sz[1]-3);

            const unsigned int k_b=ijk_lm[2*2+0]+pw;
            const unsigned int k_e=ijk_lm[2*2+1]-pw;

            const unsigned int j=(j_b-ijk_lm[2])+threadIdx.y+threadIdx.z*blockDim.y;
            const unsigned int i=(i_b-ijk_lm[0])+threadIdx.x;

            if((i>(i_e-ijk_lm[0])) || (j>(j_e-ijk_lm[2]))) return;

            const double pre_factor_6_dz = -1.0 / 64.0 / dz;

            const double smr3 = 59.0 / 48.0 * 64 * dz;
            const double smr2 = 43.0 / 48.0 * 64 * dz;
            const double smr1 = 49.0 / 48.0 * 64 * dz;
            const double spr3 = smr3;
            const double spr2 = smr2;
            const double spr1 = smr1;

            unsigned int pp=IDX_L(i, j, (k_b-ijk_lm[4]));


            for(unsigned int k=(k_b-ijk_lm[4]);k<=(k_e-ijk_lm[4]);k+=1,pp+=(tile_sz[0]*tile_sz[1]))
            {
                Du[pp] = pre_factor_6_dz *
                 (
                                        -u[pp - 3 * tile_sz[0]*tile_sz[1]]
                                        + 6.0 * u[pp - 2 * tile_sz[0]*tile_sz[1]]
                                        - 15.0 * u[pp - tile_sz[0]*tile_sz[1]]
                                        + 20.0 * u[pp]
                                        - 15.0 * u[pp + tile_sz[0]*tile_sz[1]]
                                        + 6.0 * u[pp + 2 * tile_sz[0]*tile_sz[1]]
                                        - u[pp + 3 * tile_sz[0]*tile_sz[1]]
                 );
            }

            if(ijk_lm[4]==pw)
            {


                        Du[IDX_L(i, j, (k_b-ijk_lm[4]))] = pre_factor_6_dz *
                                              (
                                                      -u[IDX_L(i, j, (k_b-ijk_lm[4]) + 4)]
                                                      + 6.0 * u[IDX_L(i, j, (k_b-ijk_lm[4]) + 3)]
                                                      - 15.0 * u[IDX_L(i, j, (k_b-ijk_lm[4]) + 2)]
                                                      + 20.0 * u[IDX_L(i, j, (k_b-ijk_lm[4]) + 1)]
                                                      - 15.0 * u[IDX_L(i, j, (k_b-ijk_lm[4]))]
                                                      + 6.0 * u[IDX_L(i, j, (k_b-ijk_lm[4]) - 1)]
                                                      - u[IDX_L(i, j, (k_b-ijk_lm[4]) - 2)]
                                              );


            }


            if(ijk_lm[5]==(sz[2]-1-pw))
            {

                        Du[IDX_L(i, j, (k_e-ijk_lm[4]))] = pre_factor_6_dz *
                                                  (
                                                          -u[IDX_L(i, j, (k_e-ijk_lm[4]) + 1)]
                                                          + 6.0 * u[IDX_L(i, j, (k_e-ijk_lm[4]))]
                                                          - 15.0 * u[IDX_L(i, j, (k_e-ijk_lm[4]) - 1)]
                                                          + 20.0 * u[IDX_L(i, j, (k_e-ijk_lm[4]) - 2)]
                                                          - 15.0 * u[IDX_L(i, j, (k_e-ijk_lm[4]) - 3)]
                                                          + 6.0 * u[IDX_L(i, j, (k_e-ijk_lm[4]) - 4)]
                                                          - u[IDX_L(i, j, (k_e-ijk_lm[4]) - 5)]
                                                  );



            }


        /*if ((bflag & (1u << OCT_DIR_BACK)) && (k_b==kb && (kz_b+ijk_lm[4])==kb) ) {

                for(unsigned int j=jy_b;j<jy_e;j++)
                    for(unsigned int i=ix_b;i<ix_e;i++)
                    {
                        Du[IDX_L(i, j, kz_b)] = (u[IDX_L(i, j, kz_b+3)]
                        - 3.0 * u[IDX_L(i, j, kz_b+2)]
                        + 3.0 * u[IDX_L(i, j, kz_b+1)]
                        - u[IDX_L(i, j, kz_b)]
                       ) / smr3;

                        Du[IDX_L(i, j, kz_b+1)] = (
                                                u[IDX_L(i, j, kz_b+4)]
                                                - 6.0 * u[IDX_L(i, j, kz_b+3)]
                                                + 12.0 * u[IDX_L(i, j, kz_b+2)]
                                                - 10.0 * u[IDX_L(i, j, kz_b+1)]
                                                + 3.0 * u[IDX_L(i, j, kz_b)]
                                        ) / smr2;
                        Du[IDX_L(i, j, kz_b+2)] = (
                                                u[IDX_L(i, j, kz_b+5)]
                                                - 6.0 * u[IDX_L(i, j, kz_b+4)]
                                                + 15.0 * u[IDX_L(i, j, kz_b+3)]
                                                - 19.0 * u[IDX_L(i, j, kz_b+2)]
                                                + 12.0 * u[IDX_L(i, j, kz_b+1)]
                                                - 3.0 * u[IDX_L(i, j, kz_b)]
                                        ) / smr1;

                    }



            }

        if ((bflag & (1u << OCT_DIR_FRONT)) && (k_e==ke && (kz_e+ijk_lm[4])==ke) ) {

                for(unsigned int j=jy_b;j<jy_e;j++)
                    for(unsigned int i=ix_b;i<ix_e;i++)
                    {

                        Du[IDX_L(i, j, kz_e - 3)] = (
                            u[IDX_L(i, j, kz_e - 6)]
                            - 6.0 * u[IDX_L(i, j, kz_e - 5)]
                            + 15.0 * u[IDX_L(i, j, kz_e - 4)]
                            - 19.0 * u[IDX_L(i, j, kz_e - 3)]
                            + 12.0 * u[IDX_L(i, j, kz_e - 2)]
                            - 3.0 * u[IDX_L(i, j, kz_e - 1)]
                    ) / spr1;

                    Du[IDX_L(i, j, kz_e - 2)] = (
                                                u[IDX_L(i, j, kz_e - 5)]
                                                - 6.0 * u[IDX_L(i, j, kz_e - 4)]
                                                + 12.0 * u[IDX_L(i, j, kz_e - 3)]
                                                - 10.0 * u[IDX_L(i, j, kz_e - 2)]
                                                + 3.0 * u[IDX_L(i, j, kz_e - 1)]
                                        ) / spr2;

                    Du[IDX_L(i, j, kz_e - 1)] = (
                                                u[IDX_L(i, j, kz_e - 4)]
                                                - 3.0 * u[IDX_L(i, j, kz_e - 3)]
                                                + 3.0 * u[IDX_L(i, j, kz_e - 2)]
                                                - u[IDX_L(i, j, kz_e - 1)]
                                        ) / spr3;

                    }



            }*/

#ifdef DEBUG_DERIVS_COMP
            if(std::isnan(Du[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
#endif


        }




} //end of namespace cuda




