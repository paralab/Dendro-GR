/**
 * @file gr with LTS support (not enabled currently)
 * @brief Main driver for cuda BSSNSolver.
 * @version 0.1
 * @date 2021-02-12
 *
 */
#include <iostream>
#include <vector>

#include "TreeNode.h"
#include "aeh.h"
#include "bssnAEH.h"
#include "bssnCtx.h"
#include "gr.h"
#include "grUtils.h"
#include "mesh.h"
#include "meshUtils.h"
#include "mpi.h"
#include "octUtils.h"
#include "parameters.h"
#include "rkBSSN.h"
#include "sdc.h"

int main(int argc, char** argv) {
    // 0- NUTS 1-UTS
    unsigned int ts_mode = 1;

    if (argc < 2) {
        std::cout << "usage: " << argv[0] << " PARAM_FILE [TS_MODE]"
                  << std::endl;
        std::cout << std::endl << "options:" << std::endl;
        std::cout << "  PARAM_FILE" << std::endl
                  << "      Path to the parameter file (.json file)"
                  << std::endl;
        std::cout << "  TS_MODE" << std::endl
                  << "      Time stepper mode." << std::endl;
        std::cout << "        0 - Spatially Adaptive Time Stepping (SATS, "
                     "Currently **NOT AVAILABLE**)"
                  << std::endl;
        std::cout << "        1 - Uniform Time Stepping (UTS, " << GRN
                  << "default" << NRM << ")" << std::endl;
        return 0;
    }

    if (argc > 2) ts_mode = std::atoi(argv[2]);

    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;

    int rank, npes;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    if (!rank) {
        std::cout << "======================================" << std::endl;
        std::cout << GRN << ":::: Now initializing BSSN Solver ::::" << NRM
                  << std::endl;
        if (ts_mode == 0) {
            std::cout << YLW
                      << "      - Running with the Non-Uniform/Spatially "
                         "Adaptive Time Stepper (NUTS/SATS)"
                      << NRM << std::endl;
        } else {
            std::cout << YLW
                      << "      - Running with the Uniform Time Stepper (UTS)"
                      << NRM << std::endl;
        }

        std::cout << "======================================" << std::endl;
    }

    // Print out CMAKE options
    if (!rank) {
#ifdef BSSN_COMPUTE_CONSTRAINTS
        std::cout << GRN << "  Compiled with BSSN_COMPUTE_CONSTRAINTS" << NRM
                  << std::endl;
#else
        std::cout << RED << "  Compiled without BSSN_COMPUTE_CONSTRAINTS" << NRM
                  << std::endl;
#endif
#ifdef BSSN_ENABLE_VTU_CONSTRAINT_OUTPUT
        std::cout << GRN << "  Compiled with BSSN_ENABLE_VTU_CONSTRAINT_OUTPUT"
                  << NRM << std::endl;
#else
        std::cout << RED
                  << "  Compiled without BSSN_ENABLE_VTU_CONSTRAINT_OUTPUT"
                  << NRM << std::endl;
#endif
#ifdef BSSN_ENABLE_VTU_OUTPUT
        std::cout << GRN << "  Compiled with BSSN_ENABLE_VTU_OUTPUT" << NRM
                  << std::endl;
#else
        std::cout << RED << "  Compiled without BSSN_ENABLE_VTU_OUTPUT" << NRM
                  << std::endl;
#endif
#ifdef BSSN_ETA_FUNCTION
        std::cout << GRN << "  Compiled with  BSSN_ETA_FUNCTION" << NRM
                  << std::endl;
#else
        std::cout << RED << "  Compiled without  BSSN_ETA_FUNCTION" << NRM
                  << std::endl;
#endif
#ifdef BSSN_EXTRACT_BH_LOCATIONS
        std::cout << GRN << "  Compiled with  BSSN_EXTRACT_BH_LOCATIONS" << NRM
                  << std::endl;
#else
        std::cout << RED << "  Compiled without  BSSN_EXTRACT_BH_LOCATIONS"
                  << NRM << std::endl;
#endif
#ifdef BSSN_EXTRACT_GRAVITATIONAL_WAVES
        std::cout << GRN << "  Compiled with  BSSN_EXTRACT_GRAVITATIONAL_WAVES"
                  << NRM << std::endl;
#else
        std::cout << RED
                  << "  Compiled without  BSSN_EXTRACT_GRAVITATIONAL_WAVES"
                  << NRM << std::endl;
#endif
#ifdef BSSN_EXTRACT_GRAVITATIONAL_WAVES
        std::cout << GRN << "  Compiled with  BSSN_EXTRACT_GRAVITATIONAL_WAVES"
                  << NRM << std::endl;
#else
        std::cout << RED
                  << "  Compiled without  BSSN_EXTRACT_GRAVITATIONAL_WAVES"
                  << NRM << std::endl;
#endif
#ifdef BSSN_GAUGE_ROCHESTER
        std::cout << GRN << "  Compiled with  BSSN_GAUGE_ROCHESTER" << NRM
                  << std::endl;
#else
        std::cout << RED << "  Compiled without  BSSN_GAUGE_ROCHESTER" << NRM
                  << std::endl;
#endif
#ifdef BSSN_KERR_SCHILD_TEST
        std::cout << GRN << "  Compiled with  BSSN_KERR_SCHILD_TEST" << NRM
                  << std::endl;
#else
        std::cout << RED << "  Compiled without  BSSN_KERR_SCHILD_TEST" << NRM
                  << std::endl;
#endif

#ifdef BSSN_REFINE_BASE_EH
        std::cout << GRN << "  Compiled with  BSSN_REFINE_BASE_EH" << NRM
                  << std::endl;
#else
        std::cout << RED << "  Compiled without  BSSN_REFINE_BASE_EH" << NRM
                  << std::endl;
#endif

#ifdef USE_FD_INTERP_FOR_UNZIP
        std::cout << GRN << "  Compiled with  USE_FD_INTERP_FOR_UNZIP" << NRM
                  << std::endl;
#else
        std::cout << RED << "  Compiled without  USE_FD_INTERP_FOR_UNZIP" << NRM
                  << std::endl;
#endif
    }

    std::vector<std::string> arg_s(argv, argv + argc);
bssn:
    printGitInformation(rank, arg_s);

    // 1 . read the parameter file.
    if (!rank) std::cout << " reading parameter file :" << argv[1] << std::endl;
    bssn::readParamFile(argv[1], comm);

    int root = std::min(1, npes - 1);
    bssn::dumpParamFile(std::cout, root, comm);

    _InitializeHcurve(bssn::BSSN_DIM);
    m_uiMaxDepth = bssn::BSSN_MAXDEPTH;

    if (bssn::BSSN_NUM_VARS % bssn::BSSN_ASYNC_COMM_K != 0) {
        if (!rank)
            std::cout << "[overlap communication error]: total BSSN_NUM_VARS: "
                      << bssn::BSSN_NUM_VARS
                      << " is not divisable by BSSN_ASYNC_COMM_K: "
                      << bssn::BSSN_ASYNC_COMM_K << std::endl;
        MPI_Abort(comm, 0);
    }

    // NOTE: this is where we originally had the check that Extract freq needed
    // to be greater than IO freq, but decoupling the parameters means that this
    // isn't "required", even though it's probably recommended.

    // 2. generate the initial grid.
    std::vector<ot::TreeNode> tmpNodes;
    std::function<void(double, double, double, double*)> f_init =
        [](double x, double y, double z, double* var) {
            bssn::initialDataFunctionWrapper(x, y, z, var);
        };
    std::function<double(double, double, double)> f_init_alpha =
        [](double x, double y, double z) {
            double var[24];
            bssn::initialDataFunctionWrapper(x, y, z, var);
            return var[0];
        };
    // std::function<void(double,double,double,double*)> f_init=[](double
    // x,double y,double z,double*var){bssn::KerrSchildData(x,y,z,var);};
    std::function<void(double, double, double, double*)> f_init_flat =
        [](double x, double y, double z, double* var) {
            bssn::minkowskiInitialData(x, y, z, var);
        };

    const unsigned int interpVars = bssn::BSSN_NUM_VARS;
    unsigned int varIndex[interpVars];
    for (unsigned int i = 0; i < bssn::BSSN_NUM_VARS; i++) varIndex[i] = i;

    /*varIndex[0]=bssn::VAR::U_ALPHA;
    varIndex[1]=bssn::VAR::U_CHI;*/
    DendroIntL localSz, globalSz;
    double t_stat;
    double t_stat_g[3];

    bssn::timer::t_f2o.start();

    if (bssn::BSSN_ENABLE_BLOCK_ADAPTIVITY) {
        if (!rank)
            std::cout << YLW << "Using block adaptive mesh. AMR disabled "
                      << NRM << std::endl;
        const Point pt_min(bssn::BSSN_BLK_MIN_X, bssn::BSSN_BLK_MIN_Y,
                           bssn::BSSN_BLK_MIN_Z);
        const Point pt_max(bssn::BSSN_BLK_MAX_X, bssn::BSSN_BLK_MAX_Y,
                           bssn::BSSN_BLK_MAX_Z);

        bssn::blockAdaptiveOctree(tmpNodes, pt_min, pt_max, m_uiMaxDepth - 2,
                                  m_uiMaxDepth, comm);
    } else {
        if (!rank)
            std::cout << YLW << "Using function2Octree. AMR enabled " << NRM
                      << std::endl;
        const unsigned int f2olmin =
            std::min(bssn::BSSN_BH1_MAX_LEV, bssn::BSSN_BH2_MAX_LEV);
        if (f2olmin < MAXDEAPTH_LEVEL_DIFF + 2) {
            if (!rank)
                std::cout << "BH min level should be larger than "
                          << (MAXDEAPTH_LEVEL_DIFF + 2) << std::endl;

            MPI_Abort(comm, 0);
        }

        std::function<void(double, double, double, double*)>* f_init_use;

        // TODO: Need to add some custom logic here to determine if
        // function2Octree should use flat initialization or
        //
        if (true) {
            f_init_use = &f_init;
        } else {
            f_init_use = &f_init_flat;
        }

        function2Octree(*f_init_use, bssn::BSSN_NUM_VARS, varIndex, interpVars,
                        tmpNodes, (f2olmin - MAXDEAPTH_LEVEL_DIFF - 2),
                        bssn::BSSN_WAVELET_TOL, bssn::BSSN_ELE_ORDER, comm);
    }

    ot::Mesh* mesh =
        ot::createMesh(tmpNodes.data(), tmpNodes.size(), bssn::BSSN_ELE_ORDER,
                       comm, 1, ot::SM_TYPE::FDM, bssn::BSSN_DENDRO_GRAIN_SZ,
                       bssn::BSSN_LOAD_IMB_TOL, bssn::BSSN_SPLIT_FIX);
    mesh->setDomainBounds(Point(bssn::BSSN_GRID_MIN_X, bssn::BSSN_GRID_MIN_Y,
                                bssn::BSSN_GRID_MIN_Z),
                          Point(bssn::BSSN_GRID_MAX_X, bssn::BSSN_GRID_MAX_Y,
                                bssn::BSSN_GRID_MAX_Z));
    unsigned int lmin, lmax;
    mesh->computeMinMaxLevel(lmin, lmax);
    if (!rank) {
        std::cout << "================= Grid Info (Before init grid "
                     "converge):==============================================="
                     "========"
                  << std::endl;
        std::cout << "lmin: " << lmin << " lmax:" << lmax << std::endl;
        std::cout << "dx: "
                  << ((bssn::BSSN_COMPD_MAX[0] - bssn::BSSN_COMPD_MIN[0]) *
                      ((1u << (m_uiMaxDepth - lmax)) /
                       ((double)bssn::BSSN_ELE_ORDER)) /
                      ((double)(1u << (m_uiMaxDepth))))
                  << std::endl;
        std::cout << "dt: "
                  << bssn::BSSN_CFL_FACTOR *
                         ((bssn::BSSN_COMPD_MAX[0] - bssn::BSSN_COMPD_MIN[0]) *
                          ((1u << (m_uiMaxDepth - lmax)) /
                           ((double)bssn::BSSN_ELE_ORDER)) /
                          ((double)(1u << (m_uiMaxDepth))))
                  << std::endl;
        std::cout << "ts mode: " << ts_mode << std::endl;
        std::cout << "========================================================="
                     "======================================================"
                  << std::endl;
    }

    // calculate the minimum dx
    bssn::BSSN_CURRENT_MIN_DX =
        ((bssn::BSSN_COMPD_MAX[0] - bssn::BSSN_COMPD_MIN[0]) *
         ((1u << (m_uiMaxDepth - lmax)) / ((double)bssn::BSSN_ELE_ORDER)) /
         ((double)(1u << (m_uiMaxDepth))));

    bssn::BSSN_RK45_TIME_STEP_SIZE =
        bssn::BSSN_CFL_FACTOR *
        ((bssn::BSSN_COMPD_MAX[0] - bssn::BSSN_COMPD_MIN[0]) *
         ((1u << (m_uiMaxDepth - lmax)) / ((double)bssn::BSSN_ELE_ORDER)) /
         ((double)(1u << (m_uiMaxDepth))));
    tmpNodes.clear();

    if (ts_mode == 1) {
        if (!rank)
            std::cout << GRN << "Now setting up the uniform time stepper!"
                      << NRM << std::endl;
        bssn::BSSNCtx* bssnCtx = new bssn::BSSNCtx(mesh);
        ts::ETS<DendroScalar, bssn::BSSNCtx>* ets =
            new ts::ETS<DendroScalar, bssn::BSSNCtx>(bssnCtx);
        ets->set_evolve_vars(bssnCtx->get_evolution_vars());

        if ((RKType)bssn::BSSN_RK_TYPE == RKType::RK3)
            ets->set_ets_coefficients(ts::ETSType::RK3);
        else if ((RKType)bssn::BSSN_RK_TYPE == RKType::RK4)
            ets->set_ets_coefficients(ts::ETSType::RK4);
        else if ((RKType)bssn::BSSN_RK_TYPE == RKType::RK45)
            ets->set_ets_coefficients(ts::ETSType::RK5);

        ets->init();
#if defined __PROFILE_CTX__ && defined __PROFILE_ETS__
        std::ofstream outfile;
        char fname[256];
        sprintf(fname, "bssnCtx_%d.txt", npes);
        if (!rank) {
            outfile.open(fname, std::ios_base::app);
            time_t now = time(0);
            // convert now to string form
            char* dt   = ctime(&now);
            outfile << "======================================================="
                       "====="
                    << std::endl;
            outfile << "Current time : " << dt << " --- " << std::endl;
            outfile << "======================================================="
                       "====="
                    << std::endl;
        }

        ets->init_pt();
        bssnCtx->reset_pt();
        ets->dump_pt(outfile);
        // bssnCtx->dump_pt(outfile);
#endif
        bool is_merge_executed = false;
        double t1              = MPI_Wtime();

        bssnaeh::initialize_aeh();

        // capture the curr step
        const DendroIntL start_step = ets->curr_step();

        while (ets->curr_time() < bssn::BSSN_RK_TIME_END) {
            const DendroIntL step            = ets->curr_step();
            const DendroScalar time          = ets->curr_time();

            bssn::BSSN_CURRENT_RK_COORD_TIME = time;
            bssn::BSSN_CURRENT_RK_STEP       = step;

            const bool isActive              = ets->is_active();
            const unsigned int rank_global   = ets->get_global_rank();

            // things that should happen **only** on time step 0
            if (step == 0) {
                if (!rank_global) {
                    std::cout << BLU
                              << "[ETS] : Timestep 0 - ensuring a few things "
                                 "are taken care of..."
                              << NRM << std::endl;
                }
                // for our scaling operation, we want to make sure that the
                // constraints are computed and handled
                bssnCtx->compute_constraint_variables();

                // make sure we write about the grid size at time 0
                bssnCtx->write_grid_summary_data();

                if (!rank_global) {
                    std::cout << BLU
                              << "[ETS] : Timestep 0 - Finished with things "
                                 "that should always be done at time 0!"
                              << NRM << std::endl;
                }
            }

            // on restore, but only on restore and not at the beginning
            if (step != 0 && step == start_step) {
                if (!rank_global)
                    std::cout << BLD << GRN
                              << "[ETS] : CHECKPOINT RESTORED. Doing "
                                 "additional cleanup.\n"
                              << NRM << std::endl;
                // there's no guarantee that the constraints will be computed
                // right on restore, due to the timing. we still need them
                // populated with "good" data
                bssnCtx->compute_constraint_variables();
            }

            const bool is_merged = bssnCtx->is_bh_merged(0.1);
            if (is_merged) {
                // make sure we set that bh is merged!
                // NOTE: don't worry, this won't update the bssn ctx object
                // after the first time, unless something modifies another
                // internal variable
                bssnCtx->set_is_merged(time, step);

                // bssn::BSSN_REMESH_TEST_FREQ=3 *
                // bssn::BSSN_REMESH_TEST_FREQ_AFTER_MERGER;
                // bssn::BSSN_MINDEPTH=5;
                // TODO: make BSSN refinement mode POST MERGER an option!

                // wkb 5 Sept 2024: disable these two lines
                // so I can test other refinement modes
                // bssn::BSSN_REFINEMENT_MODE = bssn::RefinementMode::WAMR;
                // bssn::BSSN_USE_WAVELET_TOL_FUNCTION = 1;
                bssn::BSSN_REMESH_TEST_FREQ =
                    bssn::BSSN_REMESH_TEST_FREQ_AFTER_MERGER;
                bssn::BSSN_GW_EXTRACT_FREQ =
                    bssn::BSSN_GW_EXTRACT_FREQ_AFTER_MERGER;

                // ONLY ENABLE CAKO DURING MERGER
                if (bssn::BSSN_KO_SIGMA_SCALE_BY_CONFORMAL_POST_MERGER_ONLY) {
                    bssn::BSSN_CAKO_ENABLED = true;
                }
            }

            if ((step % bssn::BSSN_REMESH_TEST_FREQ) == 0 && step != 0) {
                bool isRemesh = bssnCtx->is_remesh();
                if (isRemesh) {
                    if (!rank_global)
                        std::cout << YLW << "[ETS] : Remesh is triggered."
                                  << NRM << std::endl;

                    bssnCtx->remesh_and_gridtransfer(bssn::BSSN_DENDRO_GRAIN_SZ,
                                                     bssn::BSSN_LOAD_IMB_TOL,
                                                     bssn::BSSN_SPLIT_FIX);
                    bssn::deallocate_bssn_deriv_workspace();
                    bssn::allocate_bssn_deriv_workspace(bssnCtx->get_mesh(), 1);
                    ets->sync_with_mesh();
                    bssnCtx->calculate_full_grid_size();

                    ot::Mesh* pmesh = bssnCtx->get_mesh();
                    unsigned int lmin, lmax;
                    pmesh->computeMinMaxLevel(lmin, lmax);
                    if (!pmesh->getMPIRankGlobal())
                        printf("post merger grid level = (%d, %d)\n", lmin,
                               lmax);

                    // calculate the minimum dx
                    bssn::BSSN_CURRENT_MIN_DX =
                        ((bssn::BSSN_COMPD_MAX[0] - bssn::BSSN_COMPD_MIN[0]) *
                         ((1u << (m_uiMaxDepth - lmax)) /
                          ((double)bssn::BSSN_ELE_ORDER)) /
                         ((double)(1u << (m_uiMaxDepth))));

                    bssn::BSSN_RK45_TIME_STEP_SIZE =
                        bssn::BSSN_CFL_FACTOR *
                        ((bssn::BSSN_COMPD_MAX[0] - bssn::BSSN_COMPD_MIN[0]) *
                         ((1u << (m_uiMaxDepth - lmax)) /
                          ((double)bssn::BSSN_ELE_ORDER)) /
                         ((double)(1u << (m_uiMaxDepth))));
                    ts::TSInfo ts_in = bssnCtx->get_ts_info();
                    ts_in._m_uiTh    = bssn::BSSN_RK45_TIME_STEP_SIZE;
                    bssnCtx->set_ts_info(ts_in);

                    // REMEMBER: the true max depth of the array is two minus
                    // m_uiMaxDepth
                    if (bssn::BSSN_SCALE_VTU_AND_GW_EXTRACTION) {
                        // REMEMBER: the true max depth of the array is two
                        // minus m_uiMaxDepth
                        bssn::BSSN_IO_OUTPUT_FREQ_TRUE =
                            bssn::BSSN_IO_OUTPUT_FREQ >>
                            (m_uiMaxDepth - 2 - lmax);
                        bssn::BSSN_GW_EXTRACT_FREQ_TRUE =
                            bssn::BSSN_GW_EXTRACT_FREQ >>
                            (m_uiMaxDepth - 2 - lmax);
                        if (!rank_global)
                            std::cout << "    IO Output Freq updated to: "
                                      << bssn::BSSN_IO_OUTPUT_FREQ_TRUE
                                      << " | GW Output Freq updated to: "
                                      << bssn::BSSN_GW_EXTRACT_FREQ_TRUE
                                      << std::endl;
                    }

                    if (!rank_global) {
                        std::cout << GRN << "[ETS] : Remesh sequence finished"
                                  << NRM << std::endl;
                    }

                    // compute the constraint variables to "refresh" them on the
                    // grid for potential RHS updates
                    bssnCtx->compute_constraint_variables();
                }

                // write the grid summary data whether or not the remesh
                // happened
                bssnCtx->write_grid_summary_data();
            }

            if ((step % bssn::BSSN_TIME_STEP_OUTPUT_FREQ) == 0) {
                if (!rank_global)
                    std::cout << BLD << GRN << "[ETS - BSSN] : SOLVER UPDATE\n"
                              << NRM << "\tCurrent Step: " << ets->curr_step()
                              << "\t\tCurrent time: " << ets->curr_time()
                              << "\tdt: " << ets->ts_size() << "\t"
                              << std::endl;

                bssnCtx->terminal_output();
            }

            if ((step % bssn::BSSN_GW_EXTRACT_FREQ_TRUE) == 0) {
                if (!rank_global)
                    std::cout << "    Now extracting constraints and GW."
                              << std::endl;

                // evolving the black holes always stores the updated
                // information
                bssnCtx->evolve_bh_loc();
                bssnCtx->extract_constraints();
                bssnCtx->extract_gravitational_waves();
            }

            if ((step % bssn::BSSN_IO_OUTPUT_FREQ_TRUE) == 0) {
                // this is all IO output, except for extracting the GW waves,
                // which are "independent"

                // write to vtu, which includes writing the BH location data
                bssnCtx->write_vtu();
                bssnCtx->write_bh_coords();
            }

            if ((AEH::AEH_SOLVER_FREQ > 0) &&
                (step % AEH::AEH_SOLVER_FREQ) == 0) {
                bssnaeh::perform_aeh_step(bssnCtx, rank);
            }

            ets->evolve();

            if ((step % bssn::BSSN_CHECKPT_FREQ) == 0) {
                bssnCtx->write_checkpt();
                bssnCtx->get_mesh()->waitAll();
            }

            bssnCtx->prepare_for_next_iter();
        }

#if defined __PROFILE_CTX__ && defined __PROFILE_ETS__
        ets->dump_pt(outfile);
        // bssnCtx->dump_pt(outfile);
#endif

        double t2 = MPI_Wtime() - t1;
        double t2_g;
        par::Mpi_Allreduce(&t2, &t2_g, 1, MPI_MAX, ets->get_global_comm());
        if (!(ets->get_global_rank()))
            std::cout << " ETS time (max) : " << t2_g << std::endl;

        delete bssnCtx->get_mesh();
        delete bssnCtx;
        delete ets;

    } else {
        std::cout << RED << "Not starting solver, ts_mode needs to be set to 1!"
                  << NRM << std::endl;
    }

    MPI_Finalize();
    return 0;
}
