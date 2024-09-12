/**
 * @file nlsmInvCtx.h.h
 * @brief Inverse problems in NLSM.
 * @version 0.1
 * @date 2021-10-13
 *
 * @copyright Copyright (c) 2021
 *
 */

#pragma once
#include <iostream>

#include "daUtils.h"
#include "dendro.h"
#include "ets.h"
#include "inverseCtx.h"
#include "launcher.h"
#include "nlsmCtx.h"
#include "sdc.h"

namespace nlsm {

struct NLSM_CNTRL_VARS {
    /**@brief: Initial data Gaussian x offset */
    const double VARIATION_NLSM_ID_XC1      = 0.1;
    const double VARIATION_NLSM_ID_YC1      = 0.1;
    const double VARIATION_NLSM_ID_ZC1      = 0.1;

    /**@brief: Initial data Gaussian x offset */
    const double VARIATION_NLSM_ID_XC2      = 0.1;
    const double VARIATION_NLSM_ID_YC2      = 0.1;
    const double VARIATION_NLSM_ID_ZC2      = 0.1;

    const unsigned int NUM_WAVE_EXTRACT_LOC = 6;
    Point WAVE_EXTRAC_LOC[6]                = {Point(-1, 0, 0), Point(1, 0, 0),
                                               Point(0, -1, 0), Point(0, 1, 0),
                                               Point(0, 0, -1), Point(0, 0, 1)};

    void update_init_vars(const launcher::Launcher* la, std::ostream& sout) {
        const unsigned int split_index = la->get_comm_split_index();
        MPI_Comm local_comm            = *(la->get_sub_communicator());
        MPI_Comm global_comm           = la->get_global_communicator();

        if (split_index == 0) {
            return;

        } else if (split_index == 1) {
            nlsm::NLSM_ID_XC1 += VARIATION_NLSM_ID_XC1;

        } else if (split_index == 2) {
            nlsm::NLSM_ID_YC1 += VARIATION_NLSM_ID_YC1;

        } else if (split_index == 3) {
            nlsm::NLSM_ID_ZC1 += VARIATION_NLSM_ID_ZC1;

        } else if (split_index == 4) {
            nlsm::NLSM_ID_XC2 += VARIATION_NLSM_ID_XC2;

        } else if (split_index == 5) {
            nlsm::NLSM_ID_YC2 += VARIATION_NLSM_ID_YC2;

        } else if (split_index == 6) {
            nlsm::NLSM_ID_ZC2 += VARIATION_NLSM_ID_ZC2;
        }

        int local_rank;
        int local_npes;
        int global_rank;

        MPI_Comm_rank(local_comm, &local_rank);
        MPI_Comm_rank(global_comm, &global_rank);
        MPI_Comm_size(local_comm, &local_npes);

        if (!local_rank) {
            sout << "============= global rank : " << global_rank
                 << "=====================\n";
            sout << YLW << "\tNLSM_ID_XC1:" << nlsm::NLSM_ID_XC1 << NRM
                 << std::endl;
            sout << YLW << "\tNLSM_ID_YC1:" << nlsm::NLSM_ID_YC1 << NRM
                 << std::endl;
            sout << YLW << "\tNLSM_ID_ZC1:" << nlsm::NLSM_ID_ZC1 << NRM
                 << std::endl;
            sout << YLW << "\tNLSM_ID_XC2:" << nlsm::NLSM_ID_XC2 << NRM
                 << std::endl;
            sout << YLW << "\tNLSM_ID_YC2:" << nlsm::NLSM_ID_YC2 << NRM
                 << std::endl;
            sout << YLW << "\tNLSM_ID_ZC2:" << nlsm::NLSM_ID_ZC2 << NRM
                 << std::endl;
            sout << "=========================================================="
                    "========\n";
        }

        return;
    }
};

class InvNLSMCtx : public invp::InverseCtx<InvNLSMCtx, NLSMCtx> {
   protected:
    using invp::InverseCtx<InvNLSMCtx, NLSMCtx>::m_uiForCtx;

    using invp::InverseCtx<InvNLSMCtx, NLSMCtx>::m_uiAdjCtx;

    std::string m_uiParFile                         = "";

    ts::ETS<DendroScalar, nlsm::NLSMCtx>* m_uiTSObj = NULL;

    NLSM_CNTRL_VARS m_uiCntrlVars;

   public:
    InvNLSMCtx();

    ~InvNLSMCtx() {};

    /**@brief set parameter file name*/
    int set_parameter_filename(char* fname) {
        m_uiParFile = std::string(fname);
        return 0;
    }

    /**@brief get parameter file name*/
    const std::string& get_parameter_filename() const { return m_uiParFile; }

    /**@brief create forward problem ctx*/
    int initialize_forward_ctx();

    /**@brief create adjoint problem ctx*/
    int initialize_adjoint_ctx();

    /**@brief launch the forward problem */
    int launch_forward_solve();

    /**@brief launch the adjoint problem */
    int launch_adjoint_solve();

    /**@brief launch gradient approximation */
    int gradient_approximation();

    int extract_waves(const Point* pts, unsigned int n);
};

}  // end of namespace nlsm.
