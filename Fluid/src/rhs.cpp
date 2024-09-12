#include "rhs.h"

#include "rhshlle.h"

/*--------------------------------------------------------------------;
 *
 * RHS for relativistic fluid
 *
 *-------------------------------------------------------------------*/

// #define DEBUG_RHS

void fluid::fluid_rhs(const ot::Mesh *pMesh, DendroScalar time,
                      DendroScalar **uzipVarsRHS, const DendroScalar **uZipVars,
                      const DendroScalar **uZipVarsPrim) {
    const ot::Block *blkList   = &(*(pMesh->getLocalBlockList().begin()));
    const unsigned int numBlks = pMesh->getLocalBlockList().size();

    unsigned int offset, bflag;
    unsigned int sz[3];
    DendroScalar dx, dy, dz;
    DendroScalar ptmin[3];
    DendroScalar ptmax[3];

    const Point pt_min(fluid::FLUID_COMPD_MIN[0], fluid::FLUID_COMPD_MIN[1],
                       fluid::FLUID_COMPD_MIN[2]);
    const Point pt_max(fluid::FLUID_COMPD_MAX[0], fluid::FLUID_COMPD_MAX[1],
                       fluid::FLUID_COMPD_MAX[2]);

#ifdef DEBUG_RHS
    unsigned int totalsize = 0;
    double prim_min[5]     = {5e-10, 0.0, 0.0, 0.0, 5e-10};
    double prim_max[5]     = {0.0};
    double cons_min[5]     = {5e-10, 0.0, 0.0, 0.0, 5e-10};
    double cons_max[5]     = {0.0};
    double rhs_min[5]      = {0.0};
    double rhs_max[5]      = {0.0};
#endif
    for (unsigned int blk = 0; blk < numBlks; blk++) {
        offset = blkList[blk].getOffset();
        sz[0]  = blkList[blk].getAllocationSzX();
        sz[1]  = blkList[blk].getAllocationSzY();
        sz[2]  = blkList[blk].getAllocationSzZ();

#ifdef DEBUG_RHS
        totalsize += sz[0] * sz[1] * sz[2];
#endif

        bflag    = blkList[blk].getBlkNodeFlag();

        dx       = blkList[blk].computeDx(pt_min, pt_max);
        dy       = blkList[blk].computeDy(pt_min, pt_max);
        dz       = blkList[blk].computeDz(pt_min, pt_max);

        ptmin[0] = GRIDX_TO_X(blkList[blk].getBlockNode().minX()) - 3 * dx;
        ptmin[1] = GRIDY_TO_Y(blkList[blk].getBlockNode().minY()) - 3 * dy;
        ptmin[2] = GRIDZ_TO_Z(blkList[blk].getBlockNode().minZ()) - 3 * dz;

        ptmax[0] = GRIDX_TO_X(blkList[blk].getBlockNode().maxX()) + 3 * dx;
        ptmax[1] = GRIDY_TO_Y(blkList[blk].getBlockNode().maxY()) + 3 * dy;
        ptmax[2] = GRIDZ_TO_Z(blkList[blk].getBlockNode().maxZ()) + 3 * dz;

        // For debugging purposes, let's make sure the entire rhs block is
        // initialized to zero.
        /*for(unsigned int m = 0; m < 5; m++){
          for(unsigned int k = 0; k < sz[2]; k++){
            for(unsigned int j = 0; j < sz[1]; j++){
              for(unsigned int i = 0; i < sz[0]; i++){
                uzipVarsRHS[m][offset + i + sz[0]*(j + sz[1]*k)] = 0.0;
              }
            }
          }
        }*/

        fluid_block_rhs(time, uzipVarsRHS, uZipVars, uZipVarsPrim, offset,
                        ptmin, ptmax, sz, bflag);

#ifdef DEBUG_RHS
        // We need to make sure that the data coming out of the righthand side
        // makes sense. However, we need to exclude the ghost regions, which
        // probably just have garbage.
        for (unsigned int m = 0; m < 5; m++) {
            for (unsigned int k = 3; k < sz[2] - 3; k++) {
                for (unsigned int j = 3; j < sz[1] - 3; j++) {
                    for (unsigned int i = 3; i < sz[0] - 3; i++) {
                        unsigned int pp = i + sz[0] * (j + sz[1] * k);
                        if (uZipVars[m][offset + pp] < cons_min[m]) {
                            cons_min[m] = uZipVars[m][offset + i +
                                                      sz[0] * (j + sz[1] * k)];
                        }
                        if (uZipVars[m][offset + pp] > cons_max[m]) {
                            cons_max[m] = uZipVars[m][offset + i +
                                                      sz[0] * (j + sz[1] * k)];
                        }
                        if (uZipVarsPrim[m][offset + pp] < prim_min[m]) {
                            prim_min[m] =
                                uZipVarsPrim[m][offset + i +
                                                sz[0] * (j + sz[1] * k)];
                        }
                        if (uZipVarsPrim[m][offset + pp] > prim_max[m]) {
                            prim_max[m] =
                                uZipVarsPrim[m][offset + i +
                                                sz[0] * (j + sz[1] * k)];
                        }
                        if (uzipVarsRHS[m][offset + pp] < rhs_min[m]) {
                            rhs_min[m] =
                                uzipVarsRHS[m][offset + i +
                                               sz[0] * (j + sz[1] * k)];
                        }
                        if (uzipVarsRHS[m][offset + pp] > rhs_max[m]) {
                            rhs_max[m] =
                                uzipVarsRHS[m][offset + i +
                                               sz[0] * (j + sz[1] * k)];
                        }
                        /*if(m == 4){
                          if(uzipVarsRHS[m][offset + pp] > 100.0 ||
                        uzipVarsRHS[m][offset + pp] < -1000.0){ printf("Extreme
                        value found at (%g, %g, %g)!\n", ptmin[0] + i*dx,
                                    ptmin[1] + j*dy, ptmin[2] + k*dz);
                            printf("  rhs[4] = %g\n",uzipVarsRHS[m][offset +
                        pp]);
                          }
                        }*/
                    }
                }
            }
        }
#endif
    }

#ifdef DEBUG_RHS
    for (unsigned int i = 0; i < 5; i++) {
        // double rho_min = vecMin(uzipVarsRHS[i] +
        // blkList[0].getOffset(),totalsize); double rho_max =
        // vecMax(uzipVarsRHS[i] + blkList[0].getOffset(),totalsize);
        printf("Time t=%g, Rank %i, ||prim[%d]|| (min, max) = (%.12g, %.12g)\n",
               time, pMesh->getMPIRank(), i, prim_min[i], prim_max[i]);
        printf("Time t=%g, Rank %i, ||cons[%d]|| (min, max) = (%.12g, %.12g)\n",
               time, pMesh->getMPIRank(), i, cons_min[i], cons_max[i]);
        printf("Time t=%g, Rank %i, ||rhs[%d]|| (min, max) = (%.12g, %.12g)\n",
               time, pMesh->getMPIRank(), i, rhs_min[i], rhs_max[i]);
    }
#endif
}

void fluid::fluid_bcs_cons(const ot::Mesh *pMesh, DendroScalar **uZipVars) {
    const ot::Block *blkList   = &(*(pMesh->getLocalBlockList().begin()));
    const unsigned int numBlks = pMesh->getLocalBlockList().size();

    unsigned int offset, bflag;
    unsigned int sz[3];
    DendroScalar dx, dy, dz;
    DendroScalar ptmin[3];
    DendroScalar ptmax[3];

    const Point pt_min(fluid::FLUID_COMPD_MIN[0], fluid::FLUID_COMPD_MIN[1],
                       fluid::FLUID_COMPD_MIN[2]);
    const Point pt_max(fluid::FLUID_COMPD_MAX[0], fluid::FLUID_COMPD_MAX[1],
                       fluid::FLUID_COMPD_MAX[2]);

    for (unsigned int blk = 0; blk < numBlks; blk++) {
        offset   = blkList[blk].getOffset();
        sz[0]    = blkList[blk].getAllocationSzX();
        sz[1]    = blkList[blk].getAllocationSzY();
        sz[2]    = blkList[blk].getAllocationSzZ();

        bflag    = blkList[blk].getBlkNodeFlag();

        dx       = blkList[blk].computeDx(pt_min, pt_max);
        dy       = blkList[blk].computeDy(pt_min, pt_max);
        dz       = blkList[blk].computeDz(pt_min, pt_max);

        ptmin[0] = GRIDX_TO_X(blkList[blk].getBlockNode().minX()) - 3 * dx;
        ptmin[1] = GRIDY_TO_Y(blkList[blk].getBlockNode().minY()) - 3 * dy;
        ptmin[2] = GRIDZ_TO_Z(blkList[blk].getBlockNode().minZ()) - 3 * dz;

        ptmax[0] = GRIDX_TO_X(blkList[blk].getBlockNode().maxX()) + 3 * dx;
        ptmax[1] = GRIDY_TO_Y(blkList[blk].getBlockNode().maxY()) + 3 * dy;
        ptmax[2] = GRIDZ_TO_Z(blkList[blk].getBlockNode().maxZ()) + 3 * dz;

        fluid_block_bcs_cons(uZipVars, offset, ptmin, ptmax, sz, bflag);
    }
}

void fluid::fluid_bcs_prim(const ot::Mesh *pMesh, DendroScalar **uZipVarsPrim) {
    const ot::Block *blkList   = &(*(pMesh->getLocalBlockList().begin()));
    const unsigned int numBlks = pMesh->getLocalBlockList().size();

    unsigned int offset, bflag;
    unsigned int sz[3];
    DendroScalar dx, dy, dz;
    DendroScalar ptmin[3];
    DendroScalar ptmax[3];

    const Point pt_min(fluid::FLUID_COMPD_MIN[0], fluid::FLUID_COMPD_MIN[1],
                       fluid::FLUID_COMPD_MIN[2]);
    const Point pt_max(fluid::FLUID_COMPD_MAX[0], fluid::FLUID_COMPD_MAX[1],
                       fluid::FLUID_COMPD_MAX[2]);

    for (unsigned int blk = 0; blk < numBlks; blk++) {
        offset   = blkList[blk].getOffset();
        sz[0]    = blkList[blk].getAllocationSzX();
        sz[1]    = blkList[blk].getAllocationSzY();
        sz[2]    = blkList[blk].getAllocationSzZ();

        bflag    = blkList[blk].getBlkNodeFlag();

        dx       = blkList[blk].computeDx(pt_min, pt_max);
        dy       = blkList[blk].computeDy(pt_min, pt_max);
        dz       = blkList[blk].computeDz(pt_min, pt_max);

        ptmin[0] = GRIDX_TO_X(blkList[blk].getBlockNode().minX()) - 3 * dx;
        ptmin[1] = GRIDY_TO_Y(blkList[blk].getBlockNode().minY()) - 3 * dy;
        ptmin[2] = GRIDZ_TO_Z(blkList[blk].getBlockNode().minZ()) - 3 * dz;

        ptmax[0] = GRIDX_TO_X(blkList[blk].getBlockNode().maxX()) + 3 * dx;
        ptmax[1] = GRIDY_TO_Y(blkList[blk].getBlockNode().maxY()) + 3 * dy;
        ptmax[2] = GRIDZ_TO_Z(blkList[blk].getBlockNode().maxZ()) + 3 * dz;

        fluid_block_bcs_prim(uZipVarsPrim, offset, ptmin, ptmax, sz, bflag);
    }
}

void fluid::fluid_block_rhs(DendroScalar t, DendroScalar **unzipVarsRHS,
                            const DendroScalar **uZipVars,
                            const DendroScalar **uZipVarsPrim,
                            const unsigned int &offset,
                            const DendroScalar *pmin, const DendroScalar *pmax,
                            const unsigned int *sz, const unsigned int &bflag) {
    // Unpack all the variables.
    // Primitive variables.
    const DendroScalar *RHO  = &uZipVarsPrim[PVAR::V_RHO][offset];
    const DendroScalar *VX   = &uZipVarsPrim[PVAR::V_VX][offset];
    const DendroScalar *VY   = &uZipVarsPrim[PVAR::V_VY][offset];
    const DendroScalar *VZ   = &uZipVarsPrim[PVAR::V_VZ][offset];
    const DendroScalar *P    = &uZipVarsPrim[PVAR::V_P][offset];

    // Conserved variables.
    const DendroScalar *D    = &uZipVars[VAR::U_D][offset];
    const DendroScalar *SX   = &uZipVars[VAR::U_SX][offset];
    const DendroScalar *SY   = &uZipVars[VAR::U_SY][offset];
    const DendroScalar *SZ   = &uZipVars[VAR::U_SZ][offset];
    const DendroScalar *TAU  = &uZipVars[VAR::U_TAU][offset];

    // Right-hand side data.
    DendroScalar *D_rhs      = &unzipVarsRHS[VAR::U_D][offset];
    DendroScalar *SX_rhs     = &unzipVarsRHS[VAR::U_SX][offset];
    DendroScalar *SY_rhs     = &unzipVarsRHS[VAR::U_SY][offset];
    DendroScalar *SZ_rhs     = &unzipVarsRHS[VAR::U_SZ][offset];
    DendroScalar *TAU_rhs    = &unzipVarsRHS[VAR::U_TAU][offset];

    // Stuff these into arrays we can use in our flux and right-hand side
    // functions.
    const DendroScalar *v[5] = {RHO, VX, VY, VZ, P};
    const DendroScalar *u[5] = {D, SX, SY, SZ, TAU};
    DendroScalar *dtu[5]     = {D_rhs, SX_rhs, SY_rhs, SZ_rhs, TAU_rhs};

    const unsigned int nx    = sz[0];
    const unsigned int ny    = sz[1];
    const unsigned int nz    = sz[2];

    DendroScalar hx          = (pmax[0] - pmin[0]) / (nx - 1);
    DendroScalar hy          = (pmax[1] - pmin[1]) / (ny - 1);
    DendroScalar hz          = (pmax[2] - pmin[2]) / (nz - 1);

    int idx[3];

    unsigned int n = sz[0] * sz[1] * sz[2];

    // Normally there would be some derivatives here, but fluids are
    // obnoxious, so we vary from some of the other methods used
    // in Dendro by performing our derivatives in a method that
    // is less likely to go to pieces at the discontinuities.

    hlle::rhshlle(dtu, u, v, pmin, pmax, sz);

    // In the other codes, it looks like boundaries are applied here.
    // We delay that and apply them after the time step is performed.
    // After that, we'll convert the conserved variables to their
    // primitive counterparts by applying a root solver on a
    // transcendental equation from... well, let's not go there.
}

void fluid::fluid_block_bcs_cons(DendroScalar **uZipVars,
                                 const unsigned int &offset,
                                 const DendroScalar *pmin,
                                 const DendroScalar *pmax,
                                 const unsigned int *sz,
                                 const unsigned int &bflag) {
    const int nx = sz[0];
    const int ny = sz[1];
    const int nz = sz[2];

#pragma message("Boundary conditions (conserved vars) hardcoded for outflow.")

    DendroScalar hx = (pmax[0] - pmin[0]) / (nx - 1);
    DendroScalar hy = (pmax[1] - pmin[1]) / (ny - 1);
    DendroScalar hz = (pmax[2] - pmin[2]) / (nz - 1);

    unsigned int ib = 3;
    unsigned int jb = 3;
    unsigned int kb = 3;
    unsigned int ie = sz[0] - 3;
    unsigned int je = sz[1] - 3;
    unsigned int ke = sz[2] - 3;

    DendroScalar x, y, z;
    unsigned int pp;

    // Check if we're on the left boundary
    if (bflag & (1u << OCT_DIR_LEFT)) {
        // x = pmin[0] + ib*hx;

        for (int m = 0; m < FLUID_NUM_EVOL_VARS; m++) {
            for (unsigned int k = 0; k < nz; k++) {
                for (unsigned int j = 0; j < ny; j++) {
                    for (unsigned int i = 0; i < ib; i++) {
                        // pp = i + nx*(j + k*ny);
                        // uZipVars[m][offset+pp] = uZipVars[m][offset+ib +
                        // nx*(j + k*ny)];
                        pp = IDX(i, j, k);
                        uZipVars[m][offset + pp] =
                            uZipVars[m][offset + IDX(ib, j, k)];
                    }
                }
            }
        }
    }
    // Check if we're on the right boundary
    if (bflag & (1u << OCT_DIR_RIGHT)) {
        // x = pmin[0] + (ie-1)*hx;
        for (unsigned int m = 0; m < FLUID_NUM_EVOL_VARS; m++) {
            for (unsigned int k = 0; k < nz; k++) {
                for (unsigned int j = 0; j < ny; j++) {
                    for (unsigned int i = ie; i < nx; i++) {
                        // pp = i + nx*(j + k*ny);
                        // uZipVars[m][offset+pp] = uZipVars[m][offset+ie - 1 +
                        // nx*(j + k*ny)];
                        pp = IDX(i, j, k);
                        uZipVars[m][offset + pp] =
                            uZipVars[m][offset + IDX(ie - 1, j, k)];
                    }
                }
            }
        }
    }

    // Check if we're on the lower boundary
    if (bflag & (1u << OCT_DIR_DOWN)) {
        // y = pmin[1] + jb*hy;

        for (unsigned int m = 0; m < FLUID_NUM_EVOL_VARS; m++) {
            for (unsigned int k = 0; k < nz; k++) {
                for (unsigned int j = 0; j < jb; j++) {
                    for (unsigned int i = 0; i < nx; i++) {
                        // pp = i + nx*(j + k*ny);
                        // uZipVars[m][offset+pp] = uZipVars[m][offset+i +
                        // nx*(jb + k*ny)];
                        pp = IDX(i, j, k);
                        uZipVars[m][offset + pp] =
                            uZipVars[m][offset + IDX(i, jb, k)];
                    }
                }
            }
        }
    }
    // Check if we're on the upper boundary
    if (bflag & (1u << OCT_DIR_UP)) {
        // y = pmin[1] + (je-1)*hy;
        for (unsigned int m = 0; m < FLUID_NUM_EVOL_VARS; m++) {
            for (unsigned int k = 0; k < nz; k++) {
                for (unsigned int j = je; j < ny; j++) {
                    for (unsigned int i = 0; i < nx; i++) {
                        // pp = i + nx*(j + k*ny);
                        // uZipVars[m][offset+pp] = uZipVars[m][offset+i +
                        // nx*(je - 1 + k*ny)];
                        pp = IDX(i, j, k);
                        uZipVars[m][offset + pp] =
                            uZipVars[m][offset + IDX(i, je - 1, k)];
                    }
                }
            }
        }
    }

    // Check if we're on the back boundary
    if (bflag & (1u << OCT_DIR_BACK)) {
        // z = pmin[2] + kb*hz;
        for (unsigned int m = 0; m < FLUID_NUM_EVOL_VARS; m++) {
            for (unsigned int k = 0; k < kb; k++) {
                for (unsigned int j = 0; j < ny; j++) {
                    for (unsigned int i = 0; i < nx; i++) {
                        // pp = i + nx*(j + k*ny);
                        // uZipVars[m][offset+pp] = uZipVars[m][offset+i + nx*(j
                        // + kb*ny)];
                        pp = IDX(i, j, k);
                        uZipVars[m][offset + pp] =
                            uZipVars[m][offset + IDX(i, j, kb)];
                    }
                }
            }
        }
    }
    // Check if we're on the front boundary
    if (bflag & (1u << OCT_DIR_FRONT)) {
        // z = pmin[2] + (ke-1)*hz;
        for (unsigned int m = 0; m < FLUID_NUM_EVOL_VARS; m++) {
            for (unsigned int k = ke; k < nz; k++) {
                for (unsigned int j = 0; j < ny; j++) {
                    for (unsigned int i = 0; i < nx; i++) {
                        // pp = i + nx*(j + k*ny);
                        // uZipVars[m][offset+pp] = uZipVars[m][offset+i + nx*(j
                        // + (ke - 1)*ny)];
                        pp = IDX(i, j, k);
                        uZipVars[m][offset + pp] =
                            uZipVars[m][offset + IDX(i, j, ke - 1)];
                    }
                }
            }
        }
    }
}

void fluid::fluid_block_bcs_prim(DendroScalar **uZipVarsPrim,
                                 const unsigned int &offset,
                                 const DendroScalar *pmin,
                                 const DendroScalar *pmax,
                                 const unsigned int *sz,
                                 const unsigned int &bflag) {
    const int nx = sz[0];
    const int ny = sz[1];
    const int nz = sz[2];

#pragma message("Boundary conditions (primitive vars) hardcoded for outflow.")

    DendroScalar hx = (pmax[0] - pmin[0]) / (nx - 1);
    DendroScalar hy = (pmax[1] - pmin[1]) / (ny - 1);
    DendroScalar hz = (pmax[2] - pmin[2]) / (nz - 1);

    unsigned int ib = 3;
    unsigned int jb = 3;
    unsigned int kb = 3;
    unsigned int ie = sz[0] - 3;
    unsigned int je = sz[1] - 3;
    unsigned int ke = sz[2] - 3;

    DendroScalar x, y, z;
    unsigned int pp;

    // Check if we're on the left boundary
    if (bflag & (1u << OCT_DIR_LEFT)) {
        // x = pmin[0] + ib*hx;
        for (int m = 0; m < FLUID_NUM_PRIM_VARS; m++) {
            for (unsigned int k = 0; k < nz; k++) {
                for (unsigned int j = 0; j < ny; j++) {
                    for (unsigned int i = 0; i < ib; i++) {
                        // pp = i + nx*(j + k*ny);
                        // uZipVarsPrim[m][offset+pp] =
                        // uZipVarsPrim[m][offset+ib + nx*(j + k*ny)];
                        pp = IDX(i, j, k);
                        uZipVarsPrim[m][offset + pp] =
                            uZipVarsPrim[m][offset + IDX(ib, j, k)];
                    }
                }
            }
        }
    }
    // Check if we're on the right boundary
    if (bflag & (1u << OCT_DIR_RIGHT)) {
        // x = pmin[0] + (ie-1)*hx;

        for (unsigned int m = 0; m < FLUID_NUM_PRIM_VARS; m++) {
            for (unsigned int k = 0; k < nz; k++) {
                for (unsigned int j = 0; j < ny; j++) {
                    for (unsigned int i = ie; i < nx; i++) {
                        // pp = i + nx*(j + k*ny);
                        // uZipVarsPrim[m][offset+pp] =
                        // uZipVarsPrim[m][offset+ie - 1 + nx*(j + k*ny)];
                        pp = IDX(i, j, k);
                        uZipVarsPrim[m][offset + pp] =
                            uZipVarsPrim[m][offset + IDX(ie - 1, j, k)];
                    }
                }
            }
        }
    }

    // Check if we're on the lower boundary
    if (bflag & (1u << OCT_DIR_DOWN)) {
        // y = pmin[1] + jb*hy;

        for (unsigned int m = 0; m < FLUID_NUM_PRIM_VARS; m++) {
            for (unsigned int k = 0; k < nz; k++) {
                for (unsigned int j = 0; j < jb; j++) {
                    for (unsigned int i = 0; i < nx; i++) {
                        // pp = i + nx*(j + k*ny);
                        // uZipVarsPrim[m][offset+pp] = uZipVarsPrim[m][offset+i
                        // + nx*(jb + k*ny)];
                        pp = IDX(i, j, k);
                        uZipVarsPrim[m][offset + pp] =
                            uZipVarsPrim[m][offset + IDX(i, jb, k)];
                    }
                }
            }
        }
    }
    // Check if we're on the upper boundary
    if (bflag & (1u << OCT_DIR_UP)) {
        // y = pmin[1] + (je-1)*hy;

        for (unsigned int m = 0; m < FLUID_NUM_PRIM_VARS; m++) {
            for (unsigned int k = 0; k < nz; k++) {
                for (unsigned int j = je; j < ny; j++) {
                    for (unsigned int i = 0; i < nx; i++) {
                        // pp = i + nx*(j + k*ny);
                        // uZipVarsPrim[m][offset+pp] = uZipVarsPrim[m][offset+i
                        // + nx*(je - 1 + k*ny)];
                        pp = IDX(i, j, k);
                        uZipVarsPrim[m][offset + pp] =
                            uZipVarsPrim[m][offset + IDX(i, je - 1, k)];
                    }
                }
            }
        }
    }

    // Check if we're on the back boundary
    if (bflag & (1u << OCT_DIR_BACK)) {
        // z = pmin[2] + kb*hz;

        for (unsigned int m = 0; m < FLUID_NUM_PRIM_VARS; m++) {
            for (unsigned int k = 0; k < kb; k++) {
                for (unsigned int j = 0; j < ny; j++) {
                    for (unsigned int i = 0; i < nx; i++) {
                        // pp = i + nx*(j + k*ny);
                        // uZipVarsPrim[m][offset+pp] = uZipVarsPrim[m][offset+i
                        // + nx*(j + kb*ny)];
                        pp = IDX(i, j, k);
                        uZipVarsPrim[m][offset + pp] =
                            uZipVarsPrim[m][offset + IDX(i, j, kb)];
                    }
                }
            }
        }
    }
    // Check if we're on the front boundary
    if (bflag & (1u << OCT_DIR_FRONT)) {
        // z = pmin[2] + (ke-1)*hz;
        for (unsigned int m = 0; m < FLUID_NUM_PRIM_VARS; m++) {
            for (unsigned int k = ke; k < nz; k++) {
                for (unsigned int j = 0; j < ny; j++) {
                    for (unsigned int i = 0; i < nx; i++) {
                        // pp = i + nx*(j + k*ny);
                        // uZipVarsPrim[m][offset+pp] = uZipVarsPrim[m][offset+i
                        // + nx*(j + (ke - 1)*ny)];
                        pp = IDX(i, j, k);
                        uZipVarsPrim[m][offset + pp] =
                            uZipVarsPrim[m][offset + IDX(i, j, ke - 1)];
                    }
                }
            }
        }
    }
}
