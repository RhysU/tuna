//--------------------------------------------------------------------------
//
// Copyright (C) 2013 Rhys Ulerich
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------

#include <tuna/kernel.h>

#include "fct.h"

// Many tests rely on having a known number of outliers
enum { noutliers = 5 };

FCT_BGN()
{
    FCT_QTEST_BGN(examine)
    {
        tuna_kernel k = {};

        // Testing initial behavior before any observations provided
        fct_chk_eq_int(noutliers, sizeof(k.outliers)/sizeof(k.outliers[0]));
        fct_chk_eq_int(tuna_stats_cnt(&(k.stats)), 0);
        for (int i = 0; i < noutliers; ++i) {
            fct_chk_eq_dbl(k.outliers[i], 0);
        }

        // Record observations 1, 3, 5, 0, 7, 9, 0 checking state after each

        tuna_kernel_obs(&k, 1); // Recorded
        fct_chk_eq_int(tuna_stats_cnt(&(k.stats)), 0);
        fct_chk_eq_dbl(k.outliers[0], 0);
        fct_chk_eq_dbl(k.outliers[1], 0);
        fct_chk_eq_dbl(k.outliers[2], 0);
        fct_chk_eq_dbl(k.outliers[3], 0);
        fct_chk_eq_dbl(k.outliers[4], 1);

        tuna_kernel_obs(&k, 3); // Recorded
        fct_chk_eq_int(tuna_stats_cnt(&(k.stats)), 0);
        fct_chk_eq_dbl(k.outliers[0], 0);
        fct_chk_eq_dbl(k.outliers[1], 0);
        fct_chk_eq_dbl(k.outliers[2], 0);
        fct_chk_eq_dbl(k.outliers[3], 1);
        fct_chk_eq_dbl(k.outliers[4], 3);

        tuna_kernel_obs(&k, 5); // Recorded
        fct_chk_eq_int(tuna_stats_cnt(&(k.stats)), 0);
        fct_chk_eq_dbl(k.outliers[0], 0);
        fct_chk_eq_dbl(k.outliers[1], 0);
        fct_chk_eq_dbl(k.outliers[2], 1);
        fct_chk_eq_dbl(k.outliers[3], 3);
        fct_chk_eq_dbl(k.outliers[4], 5);

        tuna_kernel_obs(&k, 0); // Ignored
        fct_chk_eq_int(tuna_stats_cnt(&(k.stats)), 0);
        fct_chk_eq_dbl(k.outliers[0], 0);
        fct_chk_eq_dbl(k.outliers[1], 0);
        fct_chk_eq_dbl(k.outliers[2], 1);
        fct_chk_eq_dbl(k.outliers[3], 3);
        fct_chk_eq_dbl(k.outliers[4], 5);

        tuna_kernel_obs(&k, 7); // Recorded
        fct_chk_eq_int(tuna_stats_cnt(&(k.stats)), 0);
        fct_chk_eq_dbl(k.outliers[0], 0);
        fct_chk_eq_dbl(k.outliers[1], 1);
        fct_chk_eq_dbl(k.outliers[2], 3);
        fct_chk_eq_dbl(k.outliers[3], 5);
        fct_chk_eq_dbl(k.outliers[4], 7);

        tuna_kernel_obs(&k, 9); // Recorded
        fct_chk_eq_int(tuna_stats_cnt(&(k.stats)), 0);
        fct_chk_eq_dbl(k.outliers[0], 1);
        fct_chk_eq_dbl(k.outliers[1], 3);
        fct_chk_eq_dbl(k.outliers[2], 5);
        fct_chk_eq_dbl(k.outliers[3], 7);
        fct_chk_eq_dbl(k.outliers[4], 9);

        tuna_kernel_obs(&k, 0); // Ignored
        fct_chk_eq_int(tuna_stats_cnt(&(k.stats)), 0);
        fct_chk_eq_dbl(k.outliers[0], 1);
        fct_chk_eq_dbl(k.outliers[1], 3);
        fct_chk_eq_dbl(k.outliers[2], 5);
        fct_chk_eq_dbl(k.outliers[3], 7);
        fct_chk_eq_dbl(k.outliers[4], 9);
    }
    FCT_QTEST_END();

}
FCT_END()
