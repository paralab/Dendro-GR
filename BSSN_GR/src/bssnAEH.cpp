#include "bssnAEH.h"

#include "aeh.h"
#include "bssnCtx.h"

namespace bssnaeh {

std::vector<DendroScalar> hh[4];

void initialize_aeh() {
    const unsigned int ah_num_sh = (AEH::AEH_LMAX + 1) * (AEH::AEH_LMAX + 1);

    hh[0].resize(ah_num_sh, 0);
    hh[1].resize(ah_num_sh, 0);
    hh[2].resize(ah_num_sh, 0);
    hh[3].resize(ah_num_sh, 0);

    double r_plus  = 0.5 * sqrt(std::pow(TPID::par_m_plus, 2) -
                                1 * (std::pow(TPID::par_S_plus[0], 2) +
                                    std::pow(TPID::par_S_plus[1], 2) +
                                    std::pow(TPID::par_S_plus[2], 2)));
    double r_minus = 0.5 * sqrt(std::pow(TPID::par_m_minus, 2) -
                                1 * (std::pow(TPID::par_S_minus[0], 2) +
                                     std::pow(TPID::par_S_minus[1], 2) +
                                     std::pow(TPID::par_S_minus[2], 2)));

    for (unsigned int i = 0; i < hh[0].size(); i++) {
        hh[0][i] = 0.0;
        hh[1][i] = 0.0;
        hh[2][i] = 0.0;
        hh[3][i] = 0.0;
    }

    hh[0][0] = r_plus * sqrt(4 * M_PI);
    hh[2][0] = r_minus * sqrt(4 * M_PI);
}

void perform_aeh_step(bssn::BSSNCtx* const bssnCtx, const int rank) {
    const int lmax              = AEH::AEH_LMAX;
    const int ntheta            = AEH::AEH_Q_THETA;
    const int nphi              = AEH::AEH_Q_PHI;
    const unsigned int max_iter = AEH::AEH_MAXITER;

    const double rel_tol        = AEH::AEH_ATOL;
    const double abs_tol        = AEH::AEH_RTOL;
    const bool single_ah        = bssnCtx->is_bh_merged(0.1);
    const ot::Mesh* pmesh       = bssnCtx->get_mesh();

    try {
        double rlim[2] = {1e-6, 4};
        if (single_ah) {
            const double alpha = AEH::AEH_ALPHA;
            const double beta  = AEH::AEH_BETA;

            {
                // puncture-0
                const Point& bh0 = bssnCtx->get_bh0_loc();
                const Point& bh1 = bssnCtx->get_bh1_loc();

                Point bh_m_loc(0.5 * (bh0.x() + bh1.x()),
                               0.5 * (bh0.y() + bh1.y()),
                               0.5 * (bh0.z() + bh1.z()));
                aeh::SpectralAEHSolver<bssn::BSSNCtx, DendroScalar> aeh_solver(
                    bh_m_loc, bssnCtx, lmax, ntheta, nphi, rlim, false, false);
                for (unsigned int ll = 0; ll < lmax; ll += 2) {
                    aeh_solver.set_lmodes(ll);
                    aeh_solver.solve(bssnCtx, hh[2].data(), hh[3].data(),
                                     max_iter, rel_tol, abs_tol, alpha, beta,
                                     1);
                    std::swap(hh[2], hh[3]);
                }

                char fname[256];
                // output file should save in an order for finding
                // files PREFIX _
                // aeh_bh{bhid}_l{lmax}_t{ntheta}_p{nphi}
                sprintf(fname, "%s_aeh_bhMerged_l%02d_nth%03d_nphi%03d.dat",
                        bssn::BSSN_PROFILE_FILE_PREFIX.c_str(), lmax, ntheta,
                        nphi);
                aeh_solver.set_lmodes(lmax);
                aeh_solver.solve(bssnCtx, hh[2].data(), hh[3].data(), max_iter,
                                 rel_tol, abs_tol, alpha, beta, 1);
                aeh_solver.aeh_to_json(bssnCtx, hh[3].data(), fname,
                                       std::ios_base::app);
                std::swap(hh[2], hh[3]);
            }

        } else {
            {
                const double alpha = AEH::AEH_ALPHA * bssn::BH1.getBHMass();
                const double beta  = AEH::AEH_BETA;
                // puncture-0
                if (!rank) {
                    printf(
                        "------------------------------------------"
                        "------------------\n");
                    printf(
                        "------------------------ BH 1  "
                        "-----------------------------\n");
                    printf(
                        "------------------------------------------"
                        "------------------\n");
                }

                aeh::SpectralAEHSolver<bssn::BSSNCtx, DendroScalar> aeh_solver(
                    bssnCtx->get_bh0_loc(), bssnCtx, lmax, ntheta, nphi, rlim,
                    false, false);
                for (unsigned int ll = 0; ll < lmax; ll += 2) {
                    aeh_solver.set_lmodes(ll);
                    aeh_solver.solve(bssnCtx, hh[0].data(), hh[1].data(),
                                     max_iter, rel_tol, abs_tol, alpha, beta,
                                     1);
                    std::swap(hh[0], hh[1]);
                }

                char fname[256];
                sprintf(fname, "%s_aeh_bh0_l%02d_nth%03d_nphi%03d.dat",
                        bssn::BSSN_PROFILE_FILE_PREFIX.c_str(), lmax, ntheta,
                        nphi);
                aeh_solver.set_lmodes(lmax);
                aeh_solver.solve(bssnCtx, hh[0].data(), hh[1].data(), max_iter,
                                 rel_tol, abs_tol, alpha, beta, 1);
                aeh_solver.aeh_to_json(bssnCtx, hh[1].data(), fname,
                                       std::ios_base::app);
                std::swap(hh[0], hh[1]);
            }

            {
                const double alpha = AEH::AEH_ALPHA * bssn::BH2.getBHMass();
                const double beta  = AEH::AEH_BETA;
                // puncture-1
                if (!rank) {
                    printf(
                        "------------------------------------------"
                        "------------------\n");
                    printf(
                        "------------------------ BH 2  "
                        "-----------------------------\n");
                    printf(
                        "------------------------------------------"
                        "------------------\n");
                }
                aeh::SpectralAEHSolver<bssn::BSSNCtx, DendroScalar> aeh_solver(
                    bssnCtx->get_bh1_loc(), bssnCtx, lmax, ntheta, nphi, rlim,
                    false, false);
                for (unsigned int ll = 0; ll < lmax; ll += 2) {
                    aeh_solver.set_lmodes(ll);
                    aeh_solver.solve(bssnCtx, hh[2].data(), hh[3].data(),
                                     max_iter, rel_tol, abs_tol, alpha, beta,
                                     1);
                    std::swap(hh[2], hh[3]);
                }

                char fname[256];
                sprintf(fname, "%s_aeh_bh1_l%02d_nth%03d_nphi%03d.dat",
                        bssn::BSSN_PROFILE_FILE_PREFIX.c_str(), lmax, ntheta,
                        nphi);
                aeh_solver.set_lmodes(lmax);
                aeh_solver.solve(bssnCtx, hh[2].data(), hh[3].data(), max_iter,
                                 rel_tol, abs_tol, alpha, beta, 1);
                aeh_solver.aeh_to_json(bssnCtx, hh[3].data(), fname,
                                       std::ios_base::app);
                std::swap(hh[2], hh[3]);
            }
        }

        bssnCtx->get_mesh()->waitAll();
        // needed for comm growth
        par::Mpi_Bcast(hh[0].data(), hh[0].size(), 0,
                       pmesh->getMPIGlobalCommunicator());
        par::Mpi_Bcast(hh[2].data(), hh[2].size(), 0,
                       pmesh->getMPIGlobalCommunicator());
        for (unsigned int i = 0; i < hh[0].size(); i++) {
            hh[1][i] = hh[0][i];
            hh[3][i] = hh[2][i];
        }

    } catch (...) {
        bssnCtx->get_mesh()->waitAll();
        // needed for comm growth
        par::Mpi_Bcast(hh[0].data(), hh[0].size(), 0,
                       pmesh->getMPIGlobalCommunicator());
        par::Mpi_Bcast(hh[2].data(), hh[2].size(), 0,
                       pmesh->getMPIGlobalCommunicator());
        for (unsigned int i = 0; i < hh[0].size(); i++) {
            hh[1][i] = hh[0][i];
            hh[3][i] = hh[2][i];
        }
        std::cout << "Exception occurred during apparent event "
                     "horizon solver "
                  << std::endl;
    }
}

}  // namespace bssnaeh
