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
enum { noutliers = 3 };

FCT_BGN()
{
    FCT_QTEST_BGN(observe)
    {
        tuna_kernel k = {};

        // Testing initial behavior before any observations provided
        fct_chk_eq_int(noutliers, sizeof(k.outliers)/sizeof(k.outliers[0]));
        fct_chk_eq_int(tuna_stats_cnt(&(k.stats)), 0);
        for (int i = 0; i < noutliers; ++i) {
            fct_chk_eq_dbl(k.outliers[i], 0);
        }

        // Record observations 1, 3, 0, 5, 0 checking state after each
        // Until the outliers are filled with non-zeros, count stays fixed

        tuna_kernel_obs(&k, 1); // Recorded
        fct_chk_eq_int(tuna_stats_cnt(&(k.stats)), 0);
        fct_chk_eq_dbl(k.outliers[0], 0);
        fct_chk_eq_dbl(k.outliers[1], 0);
        fct_chk_eq_dbl(k.outliers[2], 1);

        tuna_kernel_obs(&k, 3); // Recorded
        fct_chk_eq_int(tuna_stats_cnt(&(k.stats)), 0);
        fct_chk_eq_dbl(k.outliers[0], 0);
        fct_chk_eq_dbl(k.outliers[1], 1);
        fct_chk_eq_dbl(k.outliers[2], 3);

        tuna_kernel_obs(&k, 0); // Ignored
        fct_chk_eq_int(tuna_stats_cnt(&(k.stats)), 0);
        fct_chk_eq_dbl(k.outliers[0], 0);
        fct_chk_eq_dbl(k.outliers[1], 1);
        fct_chk_eq_dbl(k.outliers[2], 3);

        tuna_kernel_obs(&k, 5); // Ignored
        fct_chk_eq_int(tuna_stats_cnt(&(k.stats)), 0);
        fct_chk_eq_dbl(k.outliers[0], 1);
        fct_chk_eq_dbl(k.outliers[1], 3);
        fct_chk_eq_dbl(k.outliers[2], 5);

        tuna_kernel_obs(&k, 0); // Ignored
        fct_chk_eq_int(tuna_stats_cnt(&(k.stats)), 0);
        fct_chk_eq_dbl(k.outliers[0], 1);
        fct_chk_eq_dbl(k.outliers[1], 3);
        fct_chk_eq_dbl(k.outliers[2], 5);

        // Record elapsed times that enter middle of outliers
        // and ensure that the sorted invariant remains correct

        tuna_kernel_obs(&k, 2); // Recorded 1
        fct_chk_eq_int(tuna_stats_cnt(&(k.stats)), 1);
        fct_chk_eq_dbl(tuna_stats_avg(&(k.stats)), 1);
        fct_chk_eq_dbl(k.outliers[0], 2);
        fct_chk_eq_dbl(k.outliers[1], 3);
        fct_chk_eq_dbl(k.outliers[2], 5);

        tuna_kernel_obs(&k, 6); // Recorded 2
        fct_chk_eq_int(tuna_stats_cnt(&(k.stats)),   2);
        fct_chk_eq_dbl(tuna_stats_avg(&(k.stats)), 1.5);
        fct_chk_eq_dbl(k.outliers[0], 3);
        fct_chk_eq_dbl(k.outliers[1], 5);
        fct_chk_eq_dbl(k.outliers[2], 6);

        tuna_kernel_obs(&k, 0); // Ignored
        fct_chk_eq_int(tuna_stats_cnt(&(k.stats)), 2);
        fct_chk_eq_dbl(tuna_stats_avg(&(k.stats)), 1.5);
        fct_chk_eq_dbl(k.outliers[0], 3);
        fct_chk_eq_dbl(k.outliers[1], 5);
        fct_chk_eq_dbl(k.outliers[2], 6);

        tuna_kernel_obs(&k, 2); // Recorded 2
        fct_chk_eq_int(tuna_stats_cnt(&(k.stats)), 3);
        fct_chk_eq_dbl(tuna_stats_avg(&(k.stats)), 5./3);
        fct_chk_eq_dbl(k.outliers[0], 3);
        fct_chk_eq_dbl(k.outliers[1], 5);
        fct_chk_eq_dbl(k.outliers[2], 6);

        tuna_kernel_obs(&k, 7); // Recorded 3
        fct_chk_eq_int(tuna_stats_cnt(&(k.stats)), 4);
        fct_chk_eq_dbl(tuna_stats_avg(&(k.stats)), 2);
        fct_chk_eq_dbl(k.outliers[0], 5);
        fct_chk_eq_dbl(k.outliers[1], 6);
        fct_chk_eq_dbl(k.outliers[2], 7);

    }
    FCT_QTEST_END();

}
FCT_END()
