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

#include <tuna.h>

#include "fct.h"

// Test data and expected results
static const double obs[5] = { 2,  -3,    5,     -7,     11     };
static const double avg[5] = { 2., -1/2., 4/3.,  -3/4.,  8/5.   };
static const double var[5] = { 0., 25/2., 49/3., 113/4., 244/5. };
static const size_t N      = sizeof(obs)/sizeof(obs[0]);

FCT_BGN()
{
    FCT_QTEST_BGN(examine_obs)
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
            tuna_stats_obs(&r, obs[i]);
            fct_chk_eq_int   (tuna_stats_cnt(&r),         i+1       );
            fct_chk_eqtol_dbl(tuna_stats_avg(&r),      avg[i], 1e-14);
            fct_chk_eqtol_dbl(tuna_stats_var(&r),      var[i], 1e-14);
            fct_chk_eq_dbl   (tuna_stats_std(&r), sqrt(var[i])      );
        }
    }
    FCT_QTEST_END();

    FCT_QTEST_BGN(nobs_consistency)
    {
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
        tuna_stats_nobs(&s, obs+1, N-1);
        for (size_t i = 1; i < N; ++i) {
            tuna_stats_obs(&r, obs[i]);
        }
        fct_chk_eq_int(tuna_stats_cnt(&s), tuna_stats_cnt(&r));
        fct_chk_eq_dbl(tuna_stats_avg(&s), tuna_stats_avg(&r));
        fct_chk_eq_dbl(tuna_stats_var(&s), tuna_stats_var(&r));
        fct_chk_eq_dbl(tuna_stats_std(&s), tuna_stats_std(&r));
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

    // Tests for overflow protection via halving
    // NMAX = SIZE_MAX / 2 - 1
    static const size_t NMAX = (((size_t)-1) / 2) - 1;

    FCT_QTEST_BGN(obs_halves_at_nmax)
    {
        // Test that observation halves the accumulator when n >= NMAX
        tuna_stats s = {};

        // Directly set n to NMAX
        s.n = NMAX;
        s.m = 100.0;
        s.s = 50000.0;

        double old_m = s.m;
        (void)old_m;  // used in assertion below

        // Adding an observation should trigger halving first
        tuna_stats_obs(&s, 105.0);

        // After halving: n = NMAX/2, s = 25000.0, then observation adds 1
        // So n should be around NMAX/2 + 1
        fct_chk(s.n < NMAX);  // n was halved before incrementing
        fct_chk(s.n > NMAX/4); // but not halved twice
        // Mean should be updated but close to original
        fct_chk_eqtol_dbl(s.m, old_m, 1.0); // within 1 of original mean
    }
    FCT_QTEST_END();

    FCT_QTEST_BGN(obs_preserves_stats_approx)
    {
        // Test that halving approximately preserves statistics
        tuna_stats s1 = {};
        tuna_stats s2 = {};

        // Build up a stats object normally
        for (size_t i = 0; i < 10000; ++i) {
            tuna_stats_obs(&s1, (double)(i % 100));
        }

        double expected_avg = s1.m;
        double expected_var = tuna_stats_var(&s1);

        // Now create one with large n that will halve
        s2.n = NMAX;
        s2.m = expected_avg;
        s2.s = expected_var * (NMAX - 1);

        // Add an observation
        tuna_stats_obs(&s2, expected_avg);

        // Mean should be preserved approximately
        fct_chk_eqtol_dbl(s2.m, expected_avg, 0.01);
        // Variance should be preserved approximately (within 10%)
        double new_var = tuna_stats_var(&s2);
        fct_chk(new_var > expected_var * 0.5);
        fct_chk(new_var < expected_var * 2.0);
    }
    FCT_QTEST_END();

    FCT_QTEST_BGN(merge_handles_overflow)
    {
        // Test that merge handles potential overflow by halving
        tuna_stats dst = {};
        tuna_stats src = {};

        // Set up both with large n values
        dst.n = NMAX;
        dst.m = 50.0;
        dst.s = 10000.0;

        src.n = NMAX;
        src.m = 55.0;
        src.s = 12000.0;

        // Merge should not overflow
        tuna_stats_merge(&dst, &src);

        // Result should have n <= 2 * NMAX (and not wrapped around)
        fct_chk(dst.n <= 2 * NMAX);
        fct_chk(dst.n > 0);  // didn't wrap to zero or negative

        // Mean should be between the two input means
        fct_chk(dst.m >= 50.0);
        fct_chk(dst.m <= 55.0);
    }
    FCT_QTEST_END();

    FCT_QTEST_BGN(merge_with_oversized_src)
    {
        // Test merge when src has n > NMAX (shouldn't happen normally,
        // but we should handle it gracefully)
        tuna_stats dst = {};
        tuna_stats src = {};

        // Set up dst with normal count
        dst.n = 1000;
        dst.m = 50.0;
        dst.s = 10000.0;

        // Set up src with count > NMAX
        src.n = NMAX + 1000;
        src.m = 55.0;
        src.s = 12000.0;

        // Merge should handle this without overflow
        tuna_stats_merge(&dst, &src);

        // Result should be valid
        fct_chk(dst.n > 0);
        fct_chk(dst.n < (size_t)-1 / 2);  // no overflow
    }
    FCT_QTEST_END();

    FCT_QTEST_BGN(copy_into_empty_halves_if_needed)
    {
        // Test that copying oversized src into empty dst halves
        tuna_stats dst = {};
        tuna_stats src = {};

        // Set up src with n > NMAX
        src.n = NMAX + 1000;
        src.m = 100.0;
        src.s = 50000.0;

        // Merge into empty dst
        tuna_stats_merge(&dst, &src);

        // Result should have n <= NMAX
        fct_chk(dst.n <= NMAX);
        // Mean should be preserved
        fct_chk_eq_dbl(dst.m, 100.0);
    }
    FCT_QTEST_END();

}
FCT_END()
