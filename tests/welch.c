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

#include <tuna/welch.h>

#include "fct.h"


FCT_BGN()
{
    FCT_QTEST_BGN(welch_nuinf_large1)
    {
        // Expected found using GNU Octave: welch_test (t, c, '>')
        static const double t_mu = -0.017132, t_std = 1.0047,  t_N = 10000;
        static const double c_mu =  0.094957, c_std = 0.98075, c_N =  5000;
        static const double p_expected = 1.00000;

        double p = tuna_welch1_nuinf(t_mu, t_std, t_N, c_mu, c_std, c_N);
        fct_chk_eqtol_dbl(p, p_expected, 1e-14); // FIXME Tolerance
    }
    FCT_QTEST_END();

    FCT_QTEST_BGN(welch_nuinf_large2)
    {
        // Expected found using GNU Octave: welch_test (t, c, '>')
        static const double t_mu = -0.017132,  t_std = 1.0047, t_N = 10000;
        static const double c_mu =  0.0096706, c_std = 1.0033, c_N =  5000;
        static const double p_expected = 0.93840;

        double p = tuna_welch1_nuinf(t_mu, t_std, t_N, c_mu, c_std, c_N);
        fct_chk_eqtol_dbl(p, p_expected, 1e-14); // FIXME Tolerance
    }
    FCT_QTEST_END();

    FCT_QTEST_BGN(welch_nuinf_small)
    {
        // Expected found using GNU Octave: welch_test (t, c, '>')
        static const double t_mu = 101.67, t_std = 5.6451, t_N = 6;
        static const double c_mu = 88.833, c_std = 7.1671, c_N = 6;
        static const double p_expected = 0.00339123;

        double p = tuna_welch1_nuinf(t_mu, t_std, t_N, c_mu, c_std, c_N);
        fct_chk_eqtol_dbl(p, p_expected, 1e-14); // FIXME Conservative compare
    }
    FCT_QTEST_END();
}
FCT_END()
