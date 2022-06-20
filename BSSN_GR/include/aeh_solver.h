/**
 * @brief : To solve the Apperent event horizon solver for BSSN formulation
 * @version 0.1
 * @date 2019-11-10
 * 
 * @copyright Copyright (c) 2019
 * 
 * 
 */
#pragma once

#include "mesh.h"
#include "mpi.h"
#include "grDef.h"
#include "daUtils.h"
namespace bssn
{

    enum AEHErrorType {SUCCESS, MAX_ITERATIONS_REACHED};

    /**
     * @brief :To solve the AEH using parabolic approach. 
     * @tparam T : type of the evolving var (double, float)
     * @param pMesh : underlying mesh from Dendro
     * @param s : surface normal to the Apparent Event Horizon(AEH) 3-vector (initial guess solution to $S_2$)
     * @param bssnVars : BSSN variables at a given time
     * @param tol : tolerance for the time step iterations
     * @param max_iter : maximum number of iterations
     * @return int : error code (0 for success, )
     */
    template<typename T>
    int aeh_solver(const ot::Mesh* pMesh, T** s, T** bssnVars, T tol, unsigned int max_iter)
    {
        // Notes: We might need to change the mesh during the solver to adaptively capture the mesh. (Future work)


        Point grid_limits[2];
        Point domain_limits[2];

        grid_limits[0] = Point(bssn::BSSN_OCTREE_MIN[0], bssn::BSSN_OCTREE_MIN[1], bssn::BSSN_OCTREE_MIN[2]);
        grid_limits[1] = Point(bssn::BSSN_OCTREE_MAX[0], bssn::BSSN_OCTREE_MAX[1], bssn::BSSN_OCTREE_MAX[2]);

        domain_limits[0] = Point(bssn::BSSN_COMPD_MIN[0], bssn::BSSN_COMPD_MIN[1], bssn::BSSN_COMPD_MIN[2]);
        domain_limits[1] = Point(bssn::BSSN_COMPD_MAX[0], bssn::BSSN_COMPD_MAX[1], bssn::BSSN_COMPD_MAX[2]);

        double BSSN_BH1_BOUND[] ={-1,1};
        double BSSN_BH2_BOUND[] ={-1,1};

        const unsigned int NUM_PTS_1D = 100;
        const unsigned int NUM_PTS    = NUM_PTS_1D*NUM_PTS_1D*NUM_PTS_1D;
        double GRID_RES   = (BSSN_BH1_BOUND[1]-BSSN_BH1_BOUND[0])/NUM_PTS_1D;

        std::vector<double> domain_coords_bh1;
        std::vector<double> domain_coords_bh2;

        domain_coords_bh1.resize(3*NUM_PTS);
        domain_coords_bh2.resize(3*NUM_PTS);
        
        for(unsigned int k=0; k < NUM_PTS_1D; k++)
         for(unsigned int j=0; j < NUM_PTS_1D; j++)
          for(unsigned int i=0; i < NUM_PTS_1D; i++)
          {
            unsigned int pt_id = k*NUM_PTS_1D*NUM_PTS_1D* + j*NUM_PTS_1D + i;
            domain_coords_bh1[3* (pt_id) + 0] = BSSN_BH1_BOUND[0] + i*GRID_RES;
            domain_coords_bh1[3* (pt_id) + 1] = BSSN_BH1_BOUND[0] + j*GRID_RES;
            domain_coords_bh1[3* (pt_id) + 2] = BSSN_BH1_BOUND[0] + k*GRID_RES;

            domain_coords_bh2[3* (pt_id) + 0] = BSSN_BH2_BOUND[0] + i*GRID_RES;
            domain_coords_bh2[3* (pt_id) + 1] = BSSN_BH2_BOUND[0] + j*GRID_RES;
            domain_coords_bh2[3* (pt_id) + 2] = BSSN_BH2_BOUND[0] + k*GRID_RES;

          }

        const unsigned int GATHER_ROOT_BH1=0;
        const unsigned int GATHER_ROOT_BH2=1;
        
        std::vector<T> grid_bh1_vars;
        std::vector<T> grid_bh2_vars;

        if (pMesh->isActive() && pMesh->getMPIRank() == GATHER_ROOT_BH1)
            grid_bh1_vars.resize(NUM_PTS*bssn::BSSN_NUM_VARS);

        if (pMesh->isActive() && pMesh->getMPIRank() == GATHER_ROOT_BH2)
            grid_bh2_vars.resize(NUM_PTS*bssn::BSSN_NUM_VARS);
        
        ot::da::interpolateToCoordsAndGather(pMesh,bssnVars,domain_coords_bh1.data(),NUM_PTS,grid_limits,domain_limits,grid_bh1_vars.data(),GATHER_ROOT_BH1,bssn::BSSN_NUM_VARS);
        ot::da::interpolateToCoordsAndGather(pMesh,bssnVars,domain_coords_bh2.data(),NUM_PTS,grid_limits,domain_limits,grid_bh2_vars.data(),GATHER_ROOT_BH2,bssn::BSSN_NUM_VARS);

        if (pMesh->isActive() && pMesh->getMPIRank() == GATHER_ROOT_BH1)
        {
            // launch AEH solver BH1
        }

        if (pMesh->isActive() && pMesh->getMPIRank() == GATHER_ROOT_BH2)
        {
            // launch AEH solver BH2
        }






        
    }

}

