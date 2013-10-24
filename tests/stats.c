//--------------------------------------------------------------------------
//
// Copyright (C) 2012, 2013 Rhys Ulerich
// Copyright (C) 2012, 2013 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------

#include <string.h>

#include <tuna/stats.h>

#include "fct.h"

// Test data and expected results
static const double obs[5] = { 2,  -3,    5,     -7,     11     };
static const double avg[5] = { 2., -1/2., 4/3.,  -3/4.,  8/5.   };
static const double var[5] = { 0., 25/2., 49/3., 113/4., 244/5. };
static const size_t N      = sizeof(obs)/sizeof(obs[0]);

FCT_BGN()
{
    FCT_QTEST_BGN(examine)
    {
        tuna_stats s = {};

        // Testing initial behavior before any samples provided
        fct_chk_eq_int(0U, tuna_stats_cnt(&s));
        fct_chk(isnan(tuna_stats_avg(&s)));
        fct_chk(isnan(tuna_stats_var(&s)));
        fct_chk(isnan(tuna_stats_std(&s)));

        // Test behavior after each observation
        for (size_t i = 0; i < N; ++i) {
            tuna_stats_obs(&s, obs[i]);
            fct_chk_eq_int   (tuna_stats_cnt(&s),         i+1        );
            fct_chk_eqtol_dbl(tuna_stats_avg(&s),      avg[i], 1e-14);
            fct_chk_eqtol_dbl(tuna_stats_var(&s),      var[i], 1e-14);
            fct_chk_eq_dbl   (tuna_stats_std(&s), sqrt(var[i])      );
        }

        // Clear the accumulator and retest
        memset(&s, 0, sizeof(s));

        // Testing initial behavior before any samples provided
        fct_chk_eq_int(0U, tuna_stats_cnt(&s));
        fct_chk(isnan(tuna_stats_avg(&s)));
        fct_chk(isnan(tuna_stats_var(&s)));
        fct_chk(isnan(tuna_stats_std(&s)));

        // Test behavior after each observation
        for (size_t i = 0; i < N/2; ++i) { // NB Only half used!
            tuna_stats_obs(&s, obs[i]);
            fct_chk_eq_int(tuna_stats_cnt(&s),         i+1 );
            fct_chk_eq_dbl(tuna_stats_avg(&s),      avg[i] );
            fct_chk_eq_dbl(tuna_stats_var(&s),      var[i] );
            fct_chk_eq_dbl(tuna_stats_std(&s), sqrt(var[i]));
        }

        // Copy accumulator and ensure we can resume processing
        tuna_stats r = s;

        // Test behavior after each observation
        for (size_t i = N/2; i < N; ++i) { // NB Resuming other half!
            tuna_stats_obs(&s, obs[i]);
            fct_chk_eq_int   (tuna_stats_cnt(&s),         i+1       );
            fct_chk_eqtol_dbl(tuna_stats_avg(&s),      avg[i], 1e-14);
            fct_chk_eqtol_dbl(tuna_stats_var(&s),      var[i], 1e-14);
            fct_chk_eq_dbl   (tuna_stats_std(&s), sqrt(var[i])      );
        }
    }
    FCT_QTEST_END();

    FCT_QTEST_BGN(merge_two_trivial)
    {
        // Test merge behavior on two trivial instances
        tuna_stats r1 = {};
        tuna_stats r2 = {};

        fct_chk_eq_int(0U, tuna_stats_cnt(&r1));
        fct_chk(isnan(tuna_stats_avg(&r1)));
        fct_chk(isnan(tuna_stats_var(&r1)));
        fct_chk(isnan(tuna_stats_std(&r1)));

        tuna_stats_merge(&r1, &r2);

        fct_chk_eq_int(0U, tuna_stats_cnt(&r1));
        fct_chk(isnan(tuna_stats_avg(&r1)));
        fct_chk(isnan(tuna_stats_var(&r1)));
        fct_chk(isnan(tuna_stats_std(&r1)));

    }
    FCT_QTEST_END();

    FCT_QTEST_BGN(merge_nontrivial)
    {
        tuna_stats r1, r2;

        for (size_t k = 0; k < N; ++k) {
            memset(&r1, 0, sizeof(r1));
            memset(&r2, 0, sizeof(r2));

            // Partition data, load it, merge it
            for (size_t i = 0; i < k; ++i) tuna_stats_obs(&r1, obs[i]);
            for (size_t i = k; i < N; ++i) tuna_stats_obs(&r2, obs[i]);
            fct_chk_eq_int(k,   tuna_stats_cnt(&r1));
            fct_chk_eq_int(N-k, tuna_stats_cnt(&r2));
            tuna_stats_merge(&r1, &r2);

            // Check the entire merged result matches expected
            fct_chk_eq_int   (tuna_stats_cnt(&r1),          N          );
            fct_chk_eqtol_dbl(tuna_stats_avg(&r1),      avg[N-1], 1e-14);
            fct_chk_eqtol_dbl(tuna_stats_var(&r1),      var[N-1], 1e-14);
            fct_chk_eq_dbl   (tuna_stats_std(&r1), sqrt(var[N-1])      );
        }

    }
    FCT_QTEST_END();

}
FCT_END()
