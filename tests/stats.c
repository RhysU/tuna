//--------------------------------------------------------------------------
//
// Copyright (C) 2012, 2013, 2026 Rhys Ulerich
// Copyright (C) 2012, 2013 The PECOS Development Team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------

#include <tuna.h>

#include "fct.h"

// Test data and expected results
static const double obs[5] = { 2,  -3,    5,     -7,     11     };
static const double avg[5] = { 2., -1 / 2., 4 / 3.,  -3 / 4.,  8 / 5.   };
static const double var[5] = { 0., 25 / 2., 49 / 3., 113 / 4., 244 / 5. };
static const size_t N      = sizeof(obs) / sizeof(obs[0]);

FCT_BGN()
{
    FCT_QTEST_BGN(examine_obs) {
        tuna_stats s = {};

        // Testing initial behavior before any samples provided
        fct_chk_eq_int(0U, tuna_stats_cnt(&s));
        fct_chk(isnan(tuna_stats_avg(&s)));
        fct_chk(isnan(tuna_stats_var(&s)));
        fct_chk(isnan(tuna_stats_std(&s)));

        // Test behavior after each observation
        for (size_t i = 0; i < N; ++i) {
            tuna_stats_obs(&s, obs[i]);
            fct_chk_eq_int(tuna_stats_cnt(&s),         i + 1);
            fct_chk_eqtol_dbl(tuna_stats_avg(&s),      avg[i], 1e-14);
            fct_chk_eqtol_dbl(tuna_stats_var(&s),      var[i], 1e-14);
            fct_chk_eq_dbl(tuna_stats_std(&s), sqrt(var[i]));

            // Does the behavior of tuna_stats_mom(...) match?
            tuna_stats_mom_result mom = tuna_stats_mom(&s);
            fct_chk_eq_int(tuna_stats_cnt(&s), mom.n);
            fct_chk_eq_dbl(tuna_stats_avg(&s), mom.avg);
            fct_chk_eq_dbl(tuna_stats_var(&s), mom.var);
        }

        // Clear the accumulator and retest
        memset(&s, 0, sizeof(s));

        // Testing initial behavior before any samples provided
        fct_chk_eq_int(0U, tuna_stats_cnt(&s));
        fct_chk(isnan(tuna_stats_avg(&s)));
        fct_chk(isnan(tuna_stats_var(&s)));
        fct_chk(isnan(tuna_stats_std(&s)));

        // Test behavior after each observation
        for (size_t i = 0; i < N / 2; ++i) { // NB Only half used!
            tuna_stats_obs(&s, obs[i]);
            fct_chk_eq_int(tuna_stats_cnt(&s),         i + 1);
            fct_chk_eq_dbl(tuna_stats_avg(&s),      avg[i]);
            fct_chk_eq_dbl(tuna_stats_var(&s),      var[i]);
            fct_chk_eq_dbl(tuna_stats_std(&s), sqrt(var[i]));
        }

        // Copy accumulator and ensure we can resume processing
        tuna_stats r = s;

        // Test behavior after each observation
        for (size_t i = N / 2; i < N; ++i) { // NB Resuming other half!
            tuna_stats_obs(&r, obs[i]);
            fct_chk_eq_int(tuna_stats_cnt(&r),         i + 1);
            fct_chk_eqtol_dbl(tuna_stats_avg(&r),      avg[i], 1e-14);
            fct_chk_eqtol_dbl(tuna_stats_var(&r),      var[i], 1e-14);
            fct_chk_eq_dbl(tuna_stats_std(&r), sqrt(var[i]));
        }
    }
    FCT_QTEST_END();

    FCT_QTEST_BGN(nobs_consistency) {
        tuna_stats s = {};

        // Accumulating zero observations does nothing
        tuna_stats_nobs(&s, NULL, 0);
        fct_chk_eq_int(0U, tuna_stats_cnt(&s));
        fct_chk(isnan(tuna_stats_avg(&s)));
        fct_chk(isnan(tuna_stats_var(&s)));
        fct_chk(isnan(tuna_stats_std(&s)));

        // Accumulate just one observation and check against tuna_stats_obs
        tuna_stats_nobs(&s, obs, 1);
        tuna_stats r = {};
        tuna_stats_obs(&r, obs[0]);
        fct_chk_eq_int(tuna_stats_cnt(&s), tuna_stats_cnt(&r));
        fct_chk_eq_dbl(tuna_stats_avg(&s), tuna_stats_avg(&r));
        fct_chk_eq_dbl(tuna_stats_var(&s), tuna_stats_var(&r));
        fct_chk_eq_dbl(tuna_stats_std(&s), tuna_stats_std(&r));

        // Accumulating zero observations does nothing
        tuna_stats_nobs(&s, NULL, 0);
        fct_chk_eq_int(tuna_stats_cnt(&s), tuna_stats_cnt(&r));
        fct_chk_eq_dbl(tuna_stats_avg(&s), tuna_stats_avg(&r));
        fct_chk_eq_dbl(tuna_stats_var(&s), tuna_stats_var(&r));
        fct_chk_eq_dbl(tuna_stats_std(&s), tuna_stats_std(&r));

        // Accumulate the remainder of data and check against tuna_stats_obs
        tuna_stats_nobs(&s, obs + 1, N - 1);
        for (size_t i = 1; i < N; ++i) {
            tuna_stats_obs(&r, obs[i]);
        }
        fct_chk_eq_int(tuna_stats_cnt(&s), tuna_stats_cnt(&r));
        fct_chk_eq_dbl(tuna_stats_avg(&s), tuna_stats_avg(&r));
        fct_chk_eq_dbl(tuna_stats_var(&s), tuna_stats_var(&r));
        fct_chk_eq_dbl(tuna_stats_std(&s), tuna_stats_std(&r));
    }
    FCT_QTEST_END();

    FCT_QTEST_BGN(merge_two_trivial) {
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

    FCT_QTEST_BGN(merge_nontrivial) {
        tuna_stats r1, r2;

        for (size_t k = 0; k < N; ++k) {
            memset(&r1, 0, sizeof(r1));
            memset(&r2, 0, sizeof(r2));

            // Partition data, load it, merge it
            for (size_t i = 0; i < k; ++i) {
                tuna_stats_obs(&r1, obs[i]);
            }
            for (size_t i = k; i < N; ++i) {
                tuna_stats_obs(&r2, obs[i]);
            }
            fct_chk_eq_int(k,   tuna_stats_cnt(&r1));
            fct_chk_eq_int(N - k, tuna_stats_cnt(&r2));
            tuna_stats_merge(&r1, &r2);

            // Check the entire merged result matches expected
            fct_chk_eq_int(tuna_stats_cnt(&r1),          N);
            fct_chk_eqtol_dbl(tuna_stats_avg(&r1),      avg[N - 1], 1e-14);
            fct_chk_eqtol_dbl(tuna_stats_var(&r1),      var[N - 1], 1e-14);
            fct_chk_eq_dbl(tuna_stats_std(&r1), sqrt(var[N - 1]));
        }

    }
    FCT_QTEST_END();

    FCT_QTEST_BGN(merge_self_is_nop) {
        // Merging a stats instance with itself should be a no-op
        // since merging data with itself produces no new information
        tuna_stats s = {};

        // Add some observations to create a non-trivial state
        for (size_t i = 0; i < N; ++i) {
            tuna_stats_obs(&s, obs[i]);
        }

        // Save the state before merging with itself
        const size_t cnt_before = tuna_stats_cnt(&s);
        const double avg_before = tuna_stats_avg(&s);
        const double var_before = tuna_stats_var(&s);
        const double std_before = tuna_stats_std(&s);

        // Merge the instance with itself
        tuna_stats_merge(&s, &s);

        // Verify that all statistics remain unchanged
        fct_chk_eq_int(tuna_stats_cnt(&s), cnt_before);
        fct_chk_eq_dbl(tuna_stats_avg(&s), avg_before);
        fct_chk_eq_dbl(tuna_stats_var(&s), var_before);
        fct_chk_eq_dbl(tuna_stats_std(&s), std_before);
    }
    FCT_QTEST_END();

}
FCT_END()
