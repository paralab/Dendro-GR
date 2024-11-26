#pragma once
#include <iostream>

#include "grDef.h"
#include "grUtils.h"
#include "parameters.h"
using namespace bssn;

/*----------------------------------------------------------------------;
 *
 * enforce physical constraints on BSSN variables:
 *            det(gt) = 1,  tr(At) = 0,  alpha > 0 and chi >0.
 *
 *----------------------------------------------------------------------*/
inline void enforce_bssn_constraints(double **uiVar, unsigned int node) {
    const double one_third = 1.0 / 3.0;
    double gtd[3][3], Atd[3][3];

    gtd[0][0] = uiVar[VAR::U_SYMGT0][node];
    gtd[0][1] = uiVar[VAR::U_SYMGT1][node];
    gtd[0][2] = uiVar[VAR::U_SYMGT2][node];
    gtd[1][0] = gtd[0][1];
    gtd[1][1] = uiVar[VAR::U_SYMGT3][node];
    gtd[1][2] = uiVar[VAR::U_SYMGT4][node];
    gtd[2][0] = gtd[0][2];
    gtd[2][1] = gtd[1][2];
    gtd[2][2] = uiVar[VAR::U_SYMGT5][node];

    Atd[0][0] = uiVar[VAR::U_SYMAT0][node];
    Atd[0][1] = uiVar[VAR::U_SYMAT1][node];
    Atd[0][2] = uiVar[VAR::U_SYMAT2][node];
    Atd[1][0] = Atd[0][1];
    Atd[1][1] = uiVar[VAR::U_SYMAT3][node];
    Atd[1][2] = uiVar[VAR::U_SYMAT4][node];
    Atd[2][0] = Atd[0][2];
    Atd[2][1] = Atd[1][2];
    Atd[2][2] = uiVar[VAR::U_SYMAT5][node];

    double det_gtd =
        gtd[0][0] * (gtd[1][1] * gtd[2][2] - gtd[1][2] * gtd[1][2]) -
        gtd[0][1] * gtd[0][1] * gtd[2][2] +
        2.0 * gtd[0][1] * gtd[0][2] * gtd[1][2] -
        gtd[0][2] * gtd[0][2] * gtd[1][1];

    if (det_gtd < 0.0) {
        std::cout << "metric determinent is negative " << det_gtd << std::endl;
        exit(0);
        /* FIXME: What to do here? The metric is not physical. Do we reset the
         * metric to be flat? */
        gtd[0][0] = 1.0;
        gtd[0][1] = 0.0;
        gtd[0][2] = 0.0;
        gtd[1][0] = 0.0;
        gtd[1][1] = 1.0;
        gtd[1][2] = 0.0;
        gtd[2][0] = 0.0;
        gtd[2][1] = 0.0;
        gtd[2][2] = 1.0;
        det_gtd   = 1.0;
    }
    double det_gtd_to_neg_third = 1.0 / pow(det_gtd, one_third);

    for (unsigned int j = 0; j < 3; j++) {
        for (unsigned int i = 0; i < 3; i++) {
            gtd[i][j] *= det_gtd_to_neg_third;
        }
    }

    det_gtd = gtd[0][0] * (gtd[1][1] * gtd[2][2] - gtd[1][2] * gtd[1][2]) -
              gtd[0][1] * gtd[0][1] * gtd[2][2] +
              2.0 * gtd[0][1] * gtd[0][2] * gtd[1][2] -
              gtd[0][2] * gtd[0][2] * gtd[1][1];

    double detgt_m1 = det_gtd - 1.0;

    if (fabs(detgt_m1) > 1.0e-6) {
        std::cout.precision(14);
        std::cout << "enforce_bssn_constraint: det(gtd) != 1. det="
                  << std::fixed << det_gtd << std::endl;
        std::cout << "      gtd(1,1)=" << gtd[0][0] << std::endl;
        std::cout << "      gtd(1,2)=" << gtd[0][1] << std::endl;
        std::cout << "      gtd(1,3)=" << gtd[0][2] << std::endl;
        std::cout << "      gtd(2,2)=" << gtd[1][1] << std::endl;
        std::cout << "      gtd(2,3)=" << gtd[1][2] << std::endl;
        std::cout << "      gtd(3,3)=" << gtd[2][2] << std::endl;

        exit(0);
    }

    double gtu[3][3];
    double idet_gtd = 1.0 / det_gtd;
    gtu[0][0] = idet_gtd * (gtd[1][1] * gtd[2][2] - gtd[1][2] * gtd[1][2]);
    gtu[0][1] = idet_gtd * (-gtd[0][1] * gtd[2][2] + gtd[0][2] * gtd[1][2]);
    gtu[0][2] = idet_gtd * (gtd[0][1] * gtd[1][2] - gtd[0][2] * gtd[1][1]);
    gtu[1][0] = gtu[0][1];
    gtu[1][1] = idet_gtd * (gtd[0][0] * gtd[2][2] - gtd[0][2] * gtd[0][2]);
    gtu[1][2] = idet_gtd * (-gtd[0][0] * gtd[1][2] + gtd[0][1] * gtd[0][2]);
    gtu[2][0] = gtu[0][2];
    gtu[2][1] = gtu[1][2];
    gtu[2][2] = idet_gtd * (gtd[0][0] * gtd[1][1] - gtd[0][1] * gtd[0][1]);

    /* Require Atd to be traceless. */
    double one_third_trace_Atd =
        one_third *
        (Atd[0][0] * gtu[0][0] + Atd[1][1] * gtu[1][1] + Atd[2][2] * gtu[2][2] +
         2.0 * (Atd[0][1] * gtu[0][1] + Atd[0][2] * gtu[0][2] +
                Atd[1][2] * gtu[1][2]));

    Atd[0][0] -= one_third_trace_Atd * gtd[0][0];
    Atd[0][1] -= one_third_trace_Atd * gtd[0][1];
    Atd[0][2] -= one_third_trace_Atd * gtd[0][2];
    Atd[1][1] -= one_third_trace_Atd * gtd[1][1];
    Atd[1][2] -= one_third_trace_Atd * gtd[1][2];
    Atd[2][2] -= one_third_trace_Atd * gtd[2][2];

    double tr_A = Atd[0][0] * gtu[0][0] + Atd[1][1] * gtu[1][1] +
                  Atd[2][2] * gtu[2][2] +
                  2.0 * (Atd[0][1] * gtu[0][1] + Atd[0][2] * gtu[0][2] +
                         Atd[1][2] * gtu[1][2]);

    if (fabs(tr_A) > 1.0e-6) {
        std::cout << "enforce_bssn_constraint: tr_A != 0. tr_A=" << tr_A
                  << std::endl;
        std::cout << "      Atd(1,1)=" << Atd[0][0] << std::endl;
        std::cout << "      Atd(1,2)=" << Atd[0][1] << std::endl;
        std::cout << "      Atd(1,3)=" << Atd[0][2] << std::endl;
        std::cout << "      Atd(2,2)=" << Atd[1][1] << std::endl;
        std::cout << "      Atd(2,3)=" << Atd[1][2] << std::endl;
        std::cout << "      Atd(3,3)=" << Atd[2][2] << std::endl;

        exit(0);
    }

    uiVar[VAR::U_SYMAT0][node] = Atd[0][0];
    uiVar[VAR::U_SYMAT1][node] = Atd[0][1];
    uiVar[VAR::U_SYMAT2][node] = Atd[0][2];
    uiVar[VAR::U_SYMAT3][node] = Atd[1][1];
    uiVar[VAR::U_SYMAT4][node] = Atd[1][2];
    uiVar[VAR::U_SYMAT5][node] = Atd[2][2];

    uiVar[VAR::U_SYMGT0][node] = gtd[0][0];
    uiVar[VAR::U_SYMGT1][node] = gtd[0][1];
    uiVar[VAR::U_SYMGT2][node] = gtd[0][2];
    uiVar[VAR::U_SYMGT3][node] = gtd[1][1];
    uiVar[VAR::U_SYMGT4][node] = gtd[1][2];
    uiVar[VAR::U_SYMGT5][node] = gtd[2][2];

    /* apply a floor to chi */
    if (uiVar[VAR::U_CHI][node] < CHI_FLOOR) {
        /* FIXME This needs to be fixed when we add a fluid to the code. */
        /* ! First rescale the densitized fluid variables.
            ! The include file bssn_puncture_fluid_rescale.inc
            ! must be provided in the BSSN_*MHD project.

            ! Chi must be positive to do the rescaling of fluid variables.
            if ( chi <= 0.0) {
                chi = pars.chi_floor;
            }
            else {
                // ok... go ahead and rescale the fluid variables.
            }


            */

        /* now place the floor on chi */
        uiVar[VAR::U_CHI][node] = CHI_FLOOR;
    }

    /* apply a floor to alpha */
    uiVar[VAR::U_ALPHA][node] = std::max(uiVar[VAR::U_ALPHA][node], CHI_FLOOR);
}
