//--------------------------------------------------------------------------
//
// Copyright (C) 2013 Rhys Ulerich
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------

#include <tuna.h>

#include "fct.h"

FCT_BGN()
{
    FCT_QTEST_BGN(welch1_large1)
    {
        // Expected found using GNU Octave: welch_test (t, c, '>')
        static const double t_mu = -0.017132, t_std = 1.0047,  t_N = 10000;
        static const double c_mu =  0.094957, c_std = 0.98075, c_N =  5000;
        static const double p_expected = 1.00000;

        // For sufficiently large t-distribution DOF,
        // this normal approximation should become sufficiently close.
        double p = tuna_welch1_nuinf(t_mu, t_std, t_N, c_mu, c_std, c_N);
        fct_chk_eqtol_dbl(p, p_expected, sqrt(DBL_EPSILON));

        // The real thing should do roughly as well.
        p = tuna_welch1(t_mu, t_std, t_N, c_mu, c_std, c_N);
        fct_chk_eqtol_dbl(p, p_expected, sqrt(DBL_EPSILON));
    }
    FCT_QTEST_END();

    FCT_QTEST_BGN(welch1_large2)
    {
        // Expected found using GNU Octave: welch_test (t, c, '>')
        static const double t_mu = -0.017132,  t_std = 1.0047, t_N = 10000;
        static const double c_mu =  0.0096706, c_std = 1.0033, c_N =  5000;
        static const double p_expected = 0.93840;

        // For sufficiently large t-distribution DOF,
        // this normal approximation should become sufficiently close.
        double p = tuna_welch1_nuinf(t_mu, t_std, t_N, c_mu, c_std, c_N);
        fct_chk_eqtol_dbl(p, p_expected, 5*sqrt(sqrt(DBL_EPSILON)));

        // The real thing should do roughly as well.
        p = tuna_welch1(t_mu, t_std, t_N, c_mu, c_std, c_N);
        fct_chk_eqtol_dbl(p, p_expected, 5*sqrt(sqrt(DBL_EPSILON)));
    }
    FCT_QTEST_END();

    FCT_QTEST_BGN(welch1_small)
    {
        // Expected found using GNU Octave: welch_test (t, c, '>')
        static const double t_mu = -0.35523,  t_std = 0.93662, t_N = 10;
        static const double c_mu =  0.042389, c_std = 0.87344, c_N = 10;
        static const double p_expected = 0.817989;

        // The infinite DOF estimate should be in the right ballpark.
        double p = tuna_welch1_nuinf(t_mu, t_std, t_N, c_mu, c_std, c_N);
        fct_chk_eqtol_dbl(p, p_expected, 0.01);

        // The real thing (as CDFLIB interprets it) should do better.
        p = tuna_welch1(t_mu, t_std, t_N, c_mu, c_std, c_N);
        fct_chk_eqtol_dbl(p, p_expected, 0.001);
    }
    FCT_QTEST_END();
}
FCT_END()
