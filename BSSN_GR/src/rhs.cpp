#include "rhs.h"
#include "gr.h"
#include "hadrhs.h"

using namespace std;
using namespace bssn;


void bssnRHS(double **uzipVarsRHS, const double **uZipVars, const ot::Block* blkList, unsigned int numBlocks)
{
    unsigned int offset;
    double ptmin[3], ptmax[3];
    unsigned int sz[3];
    unsigned int bflag;
    double dx,dy,dz;
    const Point pt_min(bssn::BSSN_COMPD_MIN[0],bssn::BSSN_COMPD_MIN[1],bssn::BSSN_COMPD_MIN[2]);
    const Point pt_max(bssn::BSSN_COMPD_MAX[0],bssn::BSSN_COMPD_MAX[1],bssn::BSSN_COMPD_MAX[2]);
    const unsigned int PW=bssn::BSSN_PADDING_WIDTH;


    #ifdef BSSN_ENABLE_CUDA
        cuda::BSSNComputeParams bssnParams;
        bssnParams.BSSN_LAMBDA[0]=bssn::BSSN_LAMBDA[0];
        bssnParams.BSSN_LAMBDA[1]=bssn::BSSN_LAMBDA[1];
        bssnParams.BSSN_LAMBDA[2]=bssn::BSSN_LAMBDA[2];
        bssnParams.BSSN_LAMBDA[3]=bssn::BSSN_LAMBDA[3];

        bssnParams.BSSN_LAMBDA_F[0]=bssn::BSSN_LAMBDA_F[0];
        bssnParams.BSSN_LAMBDA_F[1]=bssn::BSSN_LAMBDA_F[1];

        bssnParams.BSSN_ETA_POWER[0]=bssn::BSSN_ETA_POWER[0];
        bssnParams.BSSN_ETA_POWER[1]=bssn::BSSN_ETA_POWER[1];

        bssnParams.ETA_R0=bssn::ETA_R0;
        bssnParams.ETA_CONST=bssn::ETA_CONST;
        bssnParams.ETA_DAMPING=bssn::ETA_DAMPING;
        bssnParams.ETA_DAMPING_EXP=bssn::ETA_DAMPING_EXP;
        bssnParams.KO_DISS_SIGMA=bssn::KO_DISS_SIGMA;

        dim3 threadBlock(16,16,1);
        cuda::computeRHS(uzipVarsRHS,(const double **)uZipVars,blkList,numBlocks,(const cuda::BSSNComputeParams*) &bssnParams,threadBlock,pt_min,pt_max,1);
    #else

    for(unsigned int blk=0; blk<numBlocks; blk++)
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

        bssnrhs(uzipVarsRHS, (const double **)uZipVars, offset, ptmin, ptmax, sz, bflag);
        

    }
    #endif
}




/*----------------------------------------------------------------------;
 *
 * vector form of RHS
 *
 *----------------------------------------------------------------------*/
void bssnrhs(double **unzipVarsRHS, const double **uZipVars,
             const unsigned int& offset,
             const double *pmin, const double *pmax, const unsigned int *sz,
             const unsigned int& bflag)
{



    const double * const alpha = &uZipVars[VAR::U_ALPHA][offset];
    const double * const chi = &uZipVars[VAR::U_CHI][offset];
    const double * const K = &uZipVars[VAR::U_K][offset];
    const double * const gt0 = &uZipVars[VAR::U_SYMGT0][offset];
    const double * const gt1 = &uZipVars[VAR::U_SYMGT1][offset];
    const double * const gt2 = &uZipVars[VAR::U_SYMGT2][offset];
    const double * const gt3 = &uZipVars[VAR::U_SYMGT3][offset];
    const double * const gt4 = &uZipVars[VAR::U_SYMGT4][offset];
    const double * const gt5 = &uZipVars[VAR::U_SYMGT5][offset];
    const double * const beta0 = &uZipVars[VAR::U_BETA0][offset];
    const double * const beta1 = &uZipVars[VAR::U_BETA1][offset];
    const double * const beta2 = &uZipVars[VAR::U_BETA2][offset];
    const double * const At0 = &uZipVars[VAR::U_SYMAT0][offset];
    const double * const At1 = &uZipVars[VAR::U_SYMAT1][offset];
    const double * const At2 = &uZipVars[VAR::U_SYMAT2][offset];
    const double * const At3 = &uZipVars[VAR::U_SYMAT3][offset];
    const double * const At4 = &uZipVars[VAR::U_SYMAT4][offset];
    const double * const At5 = &uZipVars[VAR::U_SYMAT5][offset];
    const double * const Gt0 = &uZipVars[VAR::U_GT0][offset];
    const double * const Gt1 = &uZipVars[VAR::U_GT1][offset];
    const double * const Gt2 = &uZipVars[VAR::U_GT2][offset];
    const double * const B0 = &uZipVars[VAR::U_B0][offset];
    const double * const B1 = &uZipVars[VAR::U_B1][offset];
    const double * const B2 = &uZipVars[VAR::U_B2][offset];

    double * const a_rhs = &unzipVarsRHS[VAR::U_ALPHA][offset];
    double * const chi_rhs = &unzipVarsRHS[VAR::U_CHI][offset];
    double * const K_rhs = &unzipVarsRHS[VAR::U_K][offset];
    double * const gt_rhs00 = &unzipVarsRHS[VAR::U_SYMGT0][offset];
    double * const gt_rhs01 = &unzipVarsRHS[VAR::U_SYMGT1][offset];
    double * const gt_rhs02 = &unzipVarsRHS[VAR::U_SYMGT2][offset];
    double * const gt_rhs11 = &unzipVarsRHS[VAR::U_SYMGT3][offset];
    double * const gt_rhs12 = &unzipVarsRHS[VAR::U_SYMGT4][offset];
    double * const gt_rhs22 = &unzipVarsRHS[VAR::U_SYMGT5][offset];
    double * const b_rhs0 = &unzipVarsRHS[VAR::U_BETA0][offset];
    double * const b_rhs1 = &unzipVarsRHS[VAR::U_BETA1][offset];
    double * const b_rhs2 = &unzipVarsRHS[VAR::U_BETA2][offset];
    double * const At_rhs00 = &unzipVarsRHS[VAR::U_SYMAT0][offset];
    double * const At_rhs01 = &unzipVarsRHS[VAR::U_SYMAT1][offset];
    double * const At_rhs02 = &unzipVarsRHS[VAR::U_SYMAT2][offset];
    double * const At_rhs11 = &unzipVarsRHS[VAR::U_SYMAT3][offset];
    double * const At_rhs12 = &unzipVarsRHS[VAR::U_SYMAT4][offset];
    double * const At_rhs22 = &unzipVarsRHS[VAR::U_SYMAT5][offset];
    double * const Gt_rhs0 = &unzipVarsRHS[VAR::U_GT0][offset];
    double * const Gt_rhs1 = &unzipVarsRHS[VAR::U_GT1][offset];
    double * const Gt_rhs2 = &unzipVarsRHS[VAR::U_GT2][offset];
    double * const B_rhs0 = &unzipVarsRHS[VAR::U_B0][offset];
    double * const B_rhs1 = &unzipVarsRHS[VAR::U_B1][offset];
    double * const B_rhs2 = &unzipVarsRHS[VAR::U_B2][offset];

    mem::memory_pool<double>* __mem_pool = &BSSN_MEM_POOL;

    const unsigned int nx = sz[0];
    const unsigned int ny = sz[1];
    const unsigned int nz = sz[2];

    const double hx = (pmax[0] - pmin[0]) / (nx - 1);
    const double hy = (pmax[1] - pmin[1]) / (ny - 1);
    const double hz = (pmax[2] - pmin[2]) / (nz - 1);

    const unsigned int lambda[4] = {BSSN_LAMBDA[0], BSSN_LAMBDA[1],
                                    BSSN_LAMBDA[2], BSSN_LAMBDA[3]
                                   };
    const double lambda_f[2] = {BSSN_LAMBDA_F[0], BSSN_LAMBDA_F[1]};



    int idx[3];
    const unsigned int PW=bssn::BSSN_PADDING_WIDTH;
    const unsigned int n = sz[0]*sz[1]*sz[2];
    
    bssn::timer::t_deriv.start();
    
    const unsigned int BLK_SZ=n;
    double* const deriv_base = bssn::BSSN_DERIV_WORKSPACE;

// clang-format off
    #include "bssnrhs_evar_derivs.h"
    #include "bssnrhs_derivs.h"
    #include "bssnrhs_derivs_adv.h"
// clang-format on

    bssn::timer::t_deriv.stop();

    
    // loop dep. removed allowing compiler to optmize for vectorization. 
    // if (bssn::RIT_ETA_FUNCTION == 0) {
    //     // HAD eta function
    //     eta=ETA_CONST;
    //     if (r_coord >= ETA_R0) {
    //         eta *= pow( (ETA_R0/r_coord), ETA_DAMPING_EXP);
    //     }
    // }
    //cout << "begin loop" << endl;
    bssn::timer::t_rhs.start();
    for (unsigned int k = PW; k < nz-PW; k++) {
        for (unsigned int j = PW; j < ny-PW; j++) {
            #ifdef BSSN_ENABLE_AVX
                #ifdef __INTEL_COMPILER
                #pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
                #pragma ivdep
                #endif
            #endif
            for (unsigned int i = PW; i < nx-PW; i++) {

                const double x = pmin[0] + i*hx;;
                const double y = pmin[1] + j*hy;
                const double z = pmin[2] + k*hz;;
                const unsigned int pp = i + nx*(j + ny*k);
                const double r_coord = sqrt(x*x + y*y + z*z);

                const double w = r_coord / bssn::RIT_ETA_WIDTH;
                const double arg = - w*w*w*w;
                const double eta = (bssn::RIT_ETA_CENTRAL - bssn::RIT_ETA_OUTER)*exp(arg) + bssn::RIT_ETA_OUTER;

                #ifdef USE_ROCHESTER_GAUGE
                    #pragma message("BSSN: using rochester gauge")
                    #ifdef USE_ETA_FUNC
                        #pragma message("BSSN: using function eta damping")
                        #include "bssneqs_eta_func_rochester_gauge.cpp"
                    #else
                        #pragma message("BSSN: using const eta damping")
                        #include "bssneqs_eta_const_rochester_gauge.cpp"
                    #endif
                #else
                    #pragma message("BSSN: using standard gauge")
                    #ifdef USE_ETA_FUNC
                        #pragma message("BSSN: using function eta damping")
                        #include "bssneqs_eta_func_standard_gauge.cpp"
                    #else
                        #pragma message("BSSN: using const eta damping")
                        #include "bssneqs_eta_const_standard_gauge.cpp"
                    #endif

                #endif    
                
            }
        }
    }
    bssn::timer::t_rhs.stop();

    if (bflag != 0) {

        bssn::timer::t_bdyc.start();

        bssn_bcs(a_rhs, alpha, grad_0_alpha, grad_1_alpha, grad_2_alpha, pmin, pmax,
                 1.0, 1.0, sz, bflag);
        bssn_bcs(chi_rhs, chi, grad_0_chi, grad_1_chi, grad_2_chi, pmin, pmax,
                 1.0, 1.0, sz, bflag);
        bssn_bcs(K_rhs, K, grad_0_K, grad_1_K, grad_2_K, pmin, pmax,
                 1.0, 0.0, sz, bflag);

        bssn_bcs(b_rhs0, beta0, grad_0_beta0, grad_1_beta0, grad_2_beta0, pmin, pmax,
                 1.0, 0.0, sz, bflag);
        bssn_bcs(b_rhs1, beta1, grad_0_beta1, grad_1_beta1, grad_2_beta1, pmin, pmax,
                 1.0, 0.0, sz, bflag);
        bssn_bcs(b_rhs2, beta2, grad_0_beta2, grad_1_beta2, grad_2_beta2, pmin, pmax,
                 1.0, 0.0, sz, bflag);

        bssn_bcs(Gt_rhs0, Gt0, grad_0_Gt0, grad_1_Gt0, grad_2_Gt0, pmin, pmax,
                 2.0, 0.0, sz, bflag);
        bssn_bcs(Gt_rhs1, Gt1, grad_0_Gt1, grad_1_Gt1, grad_2_Gt1, pmin, pmax,
                 2.0, 0.0, sz, bflag);
        bssn_bcs(Gt_rhs2, Gt2, grad_0_Gt2, grad_1_Gt2, grad_2_Gt2, pmin, pmax,
                 2.0, 0.0, sz, bflag);

        bssn_bcs(B_rhs0, B0, grad_0_B0, grad_1_B0, grad_2_B0, pmin, pmax,
                 1.0, 0.0, sz, bflag);
        bssn_bcs(B_rhs1, B1, grad_0_B1, grad_1_B1, grad_2_B1, pmin, pmax,
                 1.0, 0.0, sz, bflag);
        bssn_bcs(B_rhs2, B2, grad_0_B2, grad_1_B2, grad_2_B2, pmin, pmax,
                 1.0, 0.0, sz, bflag);

        bssn_bcs(At_rhs00, At0, grad_0_At0, grad_1_At0, grad_2_At0, pmin, pmax,
                 2.0, 0.0, sz, bflag);
        bssn_bcs(At_rhs01, At1, grad_0_At1, grad_1_At1, grad_2_At1, pmin, pmax,
                 2.0, 0.0, sz, bflag);
        bssn_bcs(At_rhs02, At2, grad_0_At2, grad_1_At2, grad_2_At2, pmin, pmax,
                 2.0, 0.0, sz, bflag);
        bssn_bcs(At_rhs11, At3, grad_0_At3, grad_1_At3, grad_2_At3, pmin, pmax,
                 2.0, 0.0, sz, bflag);
        bssn_bcs(At_rhs12, At4, grad_0_At4, grad_1_At4, grad_2_At4, pmin, pmax,
                 2.0, 0.0, sz, bflag);
        bssn_bcs(At_rhs22, At5, grad_0_At5, grad_1_At5, grad_2_At5, pmin, pmax,
                 2.0, 0.0, sz, bflag);

        bssn_bcs(gt_rhs00, gt0, grad_0_gt0, grad_1_gt0, grad_2_gt0, pmin, pmax,
                 1.0, 1.0, sz, bflag);
        bssn_bcs(gt_rhs01, gt1, grad_0_gt1, grad_1_gt1, grad_2_gt1, pmin, pmax,
                 1.0, 0.0, sz, bflag);
        bssn_bcs(gt_rhs02, gt2, grad_0_gt2, grad_1_gt2, grad_2_gt2, pmin, pmax,
                 1.0, 0.0, sz, bflag);
        bssn_bcs(gt_rhs11, gt3, grad_0_gt3, grad_1_gt3, grad_2_gt3, pmin, pmax,
                 1.0, 1.0, sz, bflag);
        bssn_bcs(gt_rhs12, gt4, grad_0_gt4, grad_1_gt4, grad_2_gt4, pmin, pmax,
                 1.0, 0.0, sz, bflag);
        bssn_bcs(gt_rhs22, gt5, grad_0_gt5, grad_1_gt5, grad_2_gt5, pmin, pmax,
                 1.0, 1.0, sz, bflag);

        bssn::timer::t_bdyc.stop();
    }


    bssn::timer::t_deriv.start();
    #include "bssnrhs_ko_derivs.h"
    bssn::timer::t_deriv.stop();

    bssn::timer::t_rhs.start();

    const  double sigma = KO_DISS_SIGMA;


    for (unsigned int k = PW; k < nz-PW; k++) {
        for (unsigned int j = PW; j < ny-PW; j++) {
            #ifdef BSSN_ENABLE_AVX
                #ifdef __INTEL_COMPILER
                #pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
                #pragma ivdep
                #endif
            #endif
            for (unsigned int i = PW; i < nx-PW; i++) {
                const unsigned int pp = i + nx*(j + ny*k);

                a_rhs[pp]  += sigma * (grad_0_alpha[pp] + grad_1_alpha[pp] + grad_2_alpha[pp]);
                b_rhs0[pp] += sigma * (grad_0_beta0[pp] + grad_1_beta0[pp] + grad_2_beta0[pp]);
                b_rhs1[pp] += sigma * (grad_0_beta1[pp] + grad_1_beta1[pp] + grad_2_beta1[pp]);
                b_rhs2[pp] += sigma * (grad_0_beta2[pp] + grad_1_beta2[pp] + grad_2_beta2[pp]);

                gt_rhs00[pp] += sigma * (grad_0_gt0[pp] + grad_1_gt0[pp] + grad_2_gt0[pp]);
                gt_rhs01[pp] += sigma * (grad_0_gt1[pp] + grad_1_gt1[pp] + grad_2_gt1[pp]);
                gt_rhs02[pp] += sigma * (grad_0_gt2[pp] + grad_1_gt2[pp] + grad_2_gt2[pp]);
                gt_rhs11[pp] += sigma * (grad_0_gt3[pp] + grad_1_gt3[pp] + grad_2_gt3[pp]);
                gt_rhs12[pp] += sigma * (grad_0_gt4[pp] + grad_1_gt4[pp] + grad_2_gt4[pp]);
                gt_rhs22[pp] += sigma * (grad_0_gt5[pp] + grad_1_gt5[pp] + grad_2_gt5[pp]);

                chi_rhs[pp]  += sigma * (grad_0_chi[pp] + grad_1_chi[pp] + grad_2_chi[pp]);

                At_rhs00[pp] += sigma * (grad_0_At0[pp] + grad_1_At0[pp] + grad_2_At0[pp]);
                At_rhs01[pp] += sigma * (grad_0_At1[pp] + grad_1_At1[pp] + grad_2_At1[pp]);
                At_rhs02[pp] += sigma * (grad_0_At2[pp] + grad_1_At2[pp] + grad_2_At2[pp]);
                At_rhs11[pp] += sigma * (grad_0_At3[pp] + grad_1_At3[pp] + grad_2_At3[pp]);
                At_rhs12[pp] += sigma * (grad_0_At4[pp] + grad_1_At4[pp] + grad_2_At4[pp]);
                At_rhs22[pp] += sigma * (grad_0_At5[pp] + grad_1_At5[pp] + grad_2_At5[pp]);

                K_rhs[pp] += sigma * (grad_0_K[pp] + grad_1_K[pp] + grad_2_K[pp]);

                Gt_rhs0[pp] += sigma * (grad_0_Gt0[pp] + grad_1_Gt0[pp] + grad_2_Gt0[pp]);
                Gt_rhs1[pp] += sigma * (grad_0_Gt1[pp] + grad_1_Gt1[pp] + grad_2_Gt1[pp]);
                Gt_rhs2[pp] += sigma * (grad_0_Gt2[pp] + grad_1_Gt2[pp] + grad_2_Gt2[pp]);

                B_rhs0[pp] += sigma * (grad_0_B0[pp] + grad_1_B0[pp] + grad_2_B0[pp]);
                B_rhs1[pp] += sigma * (grad_0_B1[pp] + grad_1_B1[pp] + grad_2_B1[pp]);
                B_rhs2[pp] += sigma * (grad_0_B2[pp] + grad_1_B2[pp] + grad_2_B2[pp]);
            }
        }
    }

    bssn::timer::t_rhs.stop();



    bssn::timer::t_deriv.start();
    bssn::timer::t_deriv.stop();

    #if 0
        for (unsigned int m = 0; m < 24; m++) {
            std::cout<<"  || dtu("<<m<<")|| = "<<normLInfty(unzipVarsRHS[m] + offset, n)<<std::endl;
        }
    #endif



}



/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void bssn_bcs(double *f_rhs, const double *f,
              const double *dxf, const double *dyf, const double *dzf,
              const double *pmin, const double *pmax,
              const double f_falloff, const double f_asymptotic,
              const unsigned int *sz, const unsigned int &bflag)
{

    const unsigned int nx = sz[0];
    const unsigned int ny = sz[1];
    const unsigned int nz = sz[2];

    const double hx = (pmax[0] - pmin[0]) / (nx - 1);
    const double hy = (pmax[1] - pmin[1]) / (ny - 1);
    const double hz = (pmax[2] - pmin[2]) / (nz - 1);

    const unsigned int PW=bssn::BSSN_PADDING_WIDTH;

    unsigned int ib = PW;
    unsigned int jb = PW;
    unsigned int kb = PW;
    unsigned int ie = sz[0]-PW;
    unsigned int je = sz[1]-PW;
    unsigned int ke = sz[2]-PW;

    double x,y,z;
    unsigned int pp;
    double inv_r;

    //std::cout<<"boundary bssnrhs: size [ "<<nx<<", "<<ny<<", "<<nz<<" ]"<<std::endl;
    //std::cout<<"boundary bssnrhs: pmin [ "<<pmin[0]<<", "<<pmin[1]<<", "<<pmin[2]<<" ]"<<std::endl;
    //std::cout<<"boundary bssnrhs: pmax [ "<<pmax[0]<<", "<<pmax[1]<<", "<<pmax[2]<<" ]"<<std::endl;

    if (bflag & (1u<<OCT_DIR_LEFT)) {
        double x = pmin[0] + ib*hx;
        for (unsigned int k = kb; k < ke; k++) {
            z = pmin[2] + k*hz;
            for (unsigned int j = jb; j < je; j++) {
                y = pmin[1] + j*hy;
                pp = IDX(ib,j,k);
                inv_r = 1.0 / sqrt(x*x + y*y + z*z);

                f_rhs[pp] = -  inv_r * (
                                x * dxf[pp]
                                + y * dyf[pp]
                                + z * dzf[pp]
                                + f_falloff * (   f[pp] - f_asymptotic ) );

            }
        }
    }

    if (bflag & (1u<<OCT_DIR_RIGHT)) {
        x = pmin[0] + (ie-1)*hx;
        for (unsigned int k = kb; k < ke; k++) {
            z = pmin[2] + k*hz;
            for (unsigned int j = jb; j < je; j++) {
                y = pmin[1] + j*hy;
                pp = IDX((ie-1),j,k);
                inv_r = 1.0 / sqrt(x*x + y*y + z*z);

                f_rhs[pp] = -  inv_r * (
                                x * dxf[pp]
                                + y * dyf[pp]
                                + z * dzf[pp]
                                + f_falloff * (   f[pp] - f_asymptotic ) );

            }
        }
    }

    if (bflag & (1u<<OCT_DIR_DOWN)) {
        y = pmin[1] + jb*hy;
        for (unsigned int k = kb; k < ke; k++) {
            z = pmin[2] + k*hz;
            for (unsigned int i = ib; i < ie; i++) {
                x = pmin[0] + i*hx;
                inv_r = 1.0 / sqrt(x*x + y*y + z*z);
                pp = IDX(i,jb,k);

                f_rhs[pp] = -  inv_r * (
                                x * dxf[pp]
                                + y * dyf[pp]
                                + z * dzf[pp]
                                + f_falloff * (   f[pp] - f_asymptotic ) );

            }
        }
    }

    if (bflag & (1u<<OCT_DIR_UP)) {
        y = pmin[1] + (je-1)*hy;
        for (unsigned int k = kb; k < ke; k++) {
            z = pmin[2] + k*hz;
            for (unsigned int i = ib; i < ie; i++) {
                x = pmin[0] + i*hx;
                inv_r = 1.0 / sqrt(x*x + y*y + z*z);
                pp = IDX(i,(je-1),k);

                f_rhs[pp] = -  inv_r * (
                                x * dxf[pp]
                                + y * dyf[pp]
                                + z * dzf[pp]
                                + f_falloff * (   f[pp] - f_asymptotic ) );

            }
        }
    }

    if (bflag & (1u<<OCT_DIR_BACK)) {
        z = pmin[2] + kb*hz;
        for (unsigned int j = jb; j < je; j++) {
            y = pmin[1] + j*hy;
            for (unsigned int i = ib; i < ie; i++) {
                x = pmin[0] + i*hx;
                inv_r = 1.0 / sqrt(x*x + y*y + z*z);
                pp = IDX(i,j,kb);

                f_rhs[pp] = - inv_r * (
                                x * dxf[pp]
                                + y * dyf[pp]
                                + z * dzf[pp]
                                + f_falloff * (   f[pp] - f_asymptotic ) );

            }
        }
    }

    if (bflag & (1u<<OCT_DIR_FRONT)) {
        z = pmin[2] + (ke-1)*hz;
        for (unsigned int j = jb; j < je; j++) {
            y = pmin[1] + j*hy;
            for (unsigned int i = ib; i < ie; i++) {
                x = pmin[0] + i*hx;
                inv_r = 1.0 / sqrt(x*x + y*y + z*z);
                pp = IDX(i,j,(ke-1));

                f_rhs[pp] = - inv_r * (
                                x * dxf[pp]
                                + y * dyf[pp]
                                + z * dzf[pp]
                                + f_falloff * (   f[pp] - f_asymptotic ) );

            }
        }
    }

}

/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void max_spacetime_speeds( 
                           double * const lambda1max, double * const lambda2max, double * const lambda3max, 
                           const double * const alpha, 
                           const double * const beta1, const double * const beta2, const double * const beta3,
                           const double * const gtd11, const double * const gtd12, const double * const gtd13,
                           const double * const gtd22, const double * const gtd23, const double * const gtd33,
                           const double * const chi, const unsigned int *sz)
{

    const unsigned int nx = sz[0];
    const unsigned int ny = sz[1];
    const unsigned int nz = sz[2];

    const unsigned int PW=bssn::BSSN_PADDING_WIDTH;
    unsigned int ib = PW;
    unsigned int jb = PW;
    unsigned int kb = PW;
    unsigned int ie = sz[0]-PW;
    unsigned int je = sz[1]-PW;
    unsigned int ke = sz[2]-PW;

    for (unsigned int k = kb; k < ke; k++) {
        for (unsigned int j = jb; j < je; j++) {
            for (unsigned int i = ib; i < ie; i++) {
                unsigned int pp = IDX(i,j,k);
               /* note: gtu is the inverse tilde metric. It should have detgtd = 1. So, for the purposes of 
                * calculating wavespeeds, I simple set detgtd = 1. */
                double gtu11 = gtd22[pp]*gtd33[pp] - gtd23[pp]*gtd23[pp];
                double gtu22 = gtd11[pp]*gtd33[pp] - gtd13[pp]*gtd13[pp];
                double gtu33 = gtd11[pp]*gtd22[pp] - gtd12[pp]*gtd12[pp];
                if (gtu11 < 0.0 || gtu22 < 0.0 || gtu33 < 0.0) {
                    std::cout<<"Problem computing spacetime characteristics"<<std::endl;
                    std::cout<<"gtu11 = "<<gtu11<<", gtu22 = "<<gtu22<<", gtu33 = "<<gtu33<<std::endl;
                    gtu11 = 1.0; gtu22 = 1.0; gtu33 = 1.0;
                }
                double t1 = alpha[pp] * sqrt(gtu11 * chi[pp]);
                double t2 = alpha[pp] * sqrt(gtu22 * chi[pp]);
                double t3 = alpha[pp] * sqrt(gtu33 * chi[pp]);
                lambda1max[pp] = std::max( abs(-beta1[pp] + t1), abs(-beta1[pp] - t1) );
                lambda2max[pp] = std::max( abs(-beta2[pp] + t2), abs(-beta2[pp] - t2) );
                lambda3max[pp] = std::max( abs(-beta3[pp] + t3), abs(-beta3[pp] - t3) );
            }
        }
    }
 
}

/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void freeze_bcs(double *f_rhs, const unsigned int *sz, const unsigned int &bflag)
{

    const unsigned int nx = sz[0];
    const unsigned int ny = sz[1];
    const unsigned int nz = sz[2];

    const unsigned int PW=bssn::BSSN_PADDING_WIDTH;
    unsigned int ib = PW;
    unsigned int jb = PW;
    unsigned int kb = PW;
    unsigned int ie = sz[0]-PW;
    unsigned int je = sz[1]-PW;
    unsigned int ke = sz[2]-PW;

    unsigned int pp;

    if (bflag & (1u<<OCT_DIR_LEFT)) {
        for (unsigned int k = kb; k < ke; k++) {
            for (unsigned int j = jb; j < je; j++) {
                pp = IDX(ib,j,k);
                f_rhs[pp] = 0.0;
            }
        }
    }

    if (bflag & (1u<<OCT_DIR_RIGHT)) {
        for (unsigned int k = kb; k < ke; k++) {
            for (unsigned int j = jb; j < je; j++) {
                pp = IDX((ie-1),j,k);
                f_rhs[pp] = 0.0;
            }
        }
    }

    if (bflag & (1u<<OCT_DIR_DOWN)) {
        for (unsigned int k = kb; k < ke; k++) {
            for (unsigned int i = ib; i < ie; i++) {
                pp = IDX(i,jb,k);
                f_rhs[pp] = 0.0;
            }
        }
    }

    if (bflag & (1u<<OCT_DIR_UP)) {
        for (unsigned int k = kb; k < ke; k++) {
            for (unsigned int i = ib; i < ie; i++) {
                pp = IDX(i,(je-1),k);
                f_rhs[pp] = 0.0;
            }
        }
    }

    if (bflag & (1u<<OCT_DIR_BACK)) {
        for (unsigned int j = jb; j < je; j++) {
            for (unsigned int i = ib; i < ie; i++) {
                pp = IDX(i,j,kb);
                f_rhs[pp] = 0.0;
            }
        }
    }

    if (bflag & (1u<<OCT_DIR_FRONT)) {
        for (unsigned int j = jb; j < je; j++) {
            for (unsigned int i = ib; i < ie; i++) {
                pp = IDX(i,j,(ke-1));
                f_rhs[pp] = 0.0;
            }
        }
    }

}

/*----------------------------------------------------------------------;
 *
 * HAD RHS
 *
 *----------------------------------------------------------------------*/
void call_HAD_rhs()
{
    had_bssn_rhs_();
}

#if 0
/*--------------------------------------------------------------
 * Kerr-Schild data
 *--------------------------------------------------------------*/

void ks_initial_data(double x, double y, double z, double *u)
{

    u[VAR::U_ALPHA] = 0.0;
    u[VAR::U_BETA0] = 0.0;
    u[VAR::U_BETA1] = pi/5.0*cos(y)*sin(z+x);
    u[VAR::U_BETA2] = 4.0/17.0*sin(f2*x)*sin(z);

    u[VAR::U_B0] = 0.0;
    u[VAR::U_B1] = 0.0;
    u[VAR::U_B2] = 0.0;

    u[VAR::U_GT0] = Gamt_1;
    u[VAR::U_GT1] = Gamt_2;
    u[VAR::U_GT2] = Gamt_3;

    u[VAR::U_CHI] = 1.0 + exp(-4.0*cos(x)*sin(y));

    u[VAR::U_SYMGT0] = 1.00+0.2*sin(x+z)*cos(y);
    u[VAR::U_SYMGT3] = 1.00+0.2*cos(y)*cos(z+ x);
    u[VAR::U_SYMGT5] = 1.00 / ( u[VAR::U_SYMGT0] + u[VAR::U_SYMGT3]);
    u[VAR::U_SYMGT1] = 0.7*cos(x*x + y*y);
    u[VAR::U_SYMGT2] = 0.3*sin(z)*cos(x);
    u[VAR::U_SYMGT4] = -0.5*sin(x*x)*cos(y)*cos(z);

    u[VAR::U_K] = 5.0*exp(-4.0*cos(x)*sin(y))/(5.0+sin(x))*cos(x)
                  +5.0*exp(-4.0*cos(x)*sin(y))/(5.0+cos(y))*cos(y)
                  +0.4*(25.0+5.0*cos(y)+5.0*sin(x)+sin(x)*cos(y))
                  *exp(-4.0*cos(x)*sin(y))*cos(z);

    u[VAR::U_SYMAT0] = exp(-4.0*cos(x)*sin(y))*(cos(x) -0.3333333333*exp(4.0*cos(x)*sin(y)) *(1.0+0.2*sin(x))*(5.0*exp(-4.0*cos(x)*sin(y)) /(5.0+sin(x))*cos(x)+5.0*exp(-4.0*cos(x)*sin(y)) /(5.0+cos(y))*cos(y)+0.04*(25.0+5.0*cos(y) +5.0*sin(x)+sin(x)*cos(y))*exp(-4.0*cos(x)*sin(y))*cos(z)));
    u[VAR::U_SYMAT1] = 1.0 + x*z/(0.1 + x*x + y*y + z*z);
    u[VAR::U_SYMAT2] = 1.3 - x*y/(3.0 + x*x + 2.0*y*y + z*z)*(x*x+z*z);
    u[VAR::U_SYMAT3] = exp(-4.0*cos(x)*sin(y))*(cos(y)-0.33333333330*exp(4*cos(x)*sin(y))*(1+0.2*cos(y))*(5.0*exp(-4.0*cos(x)*sin(y))/(5.0+sin(x))*cos(x)+5.0*exp(-4.0*cos(x)*sin(y))/(5.0+cos(y))*cos(y)+0.04*(25.0+5.0*cos(y)+5.0*sin(x)+sin(x)*cos(y))*exp(-4.0*cos(x)*sin(y))*cos(z)));
    u[VAR::U_SYMAT4] = -1.0 + y*z/(1.0 + 3.0*x*x + y*y + z*z);
    u[VAR::U_SYMAT5] = exp(-4.0*cos(x)*sin(y))*(cos(z)-0.3333333333*exp(4*cos(x)*sin(y))/(1+0.2*sin(x))/(1+0.2*cos(y))*(5.0*exp(-4.0*cos(x)*sin(y))/(5.0+sin(x))*cos(x)+5.0*exp(-4.0*cos(x)*sin(y))/(5.0+cos(y))*cos(y)+0.04*(25.0+5.0*cos(y)+5.0*sin(x)+sin(x)*cos(y))*exp(-4.0*cos(x)*sin(y))*cos(z)));

}
#endif
/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
