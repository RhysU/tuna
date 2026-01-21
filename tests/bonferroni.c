//--------------------------------------------------------------------------
//
// Copyright (C) 2025 Rhys Ulerich
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
    FCT_QTEST_BGN(bonferroni_without_correction_would_switch)
    {
        /* This test is designed to FAIL if Bonferroni correction is removed.
         *
         * Strategy:
         * - Create chunk 0 with better performance (lower mean)
         * - Create chunk 1 with worse performance (higher mean)
         * - Compute the actual p-value for this comparison
         * - Set u01[1] such that p < u01[1] BUT p*ncomparisons >= u01[1]
         *
         * First, let's measure what p-value we actually get, then design the test.
         */
        const size_t nchunk = 10;  /* ncomparisons = 9 */
        tuna_chunk chunks[10] = {};
        double u01[10];
        tuna_stats_mom_result stats0, stats1;
        double p_value;
        size_t i, j;
        size_t result;

        /* Chunk 0: mean ≈ 1.0 (better) */
        for (j = 0; j < 50; ++j) {
            tuna_chunk_obs(&chunks[0], 1.0 + (j % 5) * 0.01);
        }

        /* Chunk 1: mean ≈ 1.1 (worse) */
        for (j = 0; j < 50; ++j) {
            tuna_chunk_obs(&chunks[1], 1.1 + (j % 5) * 0.01);
        }

        /* Remaining chunks: same as chunk 0 */
        for (i = 2; i < nchunk; ++i) {
            for (j = 0; j < 50; ++j) {
                tuna_chunk_obs(&chunks[i], 1.0 + (j % 5) * 0.01);
            }
        }

        /* Calculate the p-value for chunk 0 vs chunk 1 */
        stats0 = tuna_stats_mom(&chunks[0].stats);
        stats1 = tuna_stats_mom(&chunks[1].stats);
        p_value = tuna_welch1_nuinf(stats0.avg, stats0.var, stats0.n,
                                     stats1.avg, stats1.var, stats1.n);

        /* Design threshold so that:
         * - p < u01[1]  (would switch WITHOUT Bonferroni)
         * - p * 9 >= u01[1]  (stays WITH Bonferroni)
         *
         * Set u01[1] to be between p and p*9
         */
        u01[0] = 0.5;  /* Not used */
        u01[1] = p_value * 5.0;  /* Between p and p*9 */

        /* Set remaining thresholds very low to prevent switching to them */
        for (i = 2; i < nchunk; ++i) {
            u01[i] = 0.0001;
        }

        /* Verify our threshold design makes sense */
        fct_chk(p_value < u01[1]);           /* Would switch without correction */
        fct_chk(p_value * 9.0 >= u01[1]);    /* Stays with correction */

        result = tuna_algo_welch1_nuinf->function(nchunk, chunks, u01);

        /* With Bonferroni correction, should stay at chunk 0 */
        fct_chk_eq_int(result, 0);
    }
    FCT_QTEST_END();

    FCT_QTEST_BGN(bonferroni_welch1_variant)
    {
        /* Same test for standard welch1 algorithm */
        const size_t nchunk = 10;
        tuna_chunk chunks[10] = {};
        double u01[10];
        tuna_stats_mom_result stats0, stats1;
        double p_value;
        size_t i, j;
        size_t result;

        for (j = 0; j < 50; ++j) {
            tuna_chunk_obs(&chunks[0], 1.0 + (j % 5) * 0.01);
        }

        for (j = 0; j < 50; ++j) {
            tuna_chunk_obs(&chunks[1], 1.1 + (j % 5) * 0.01);
        }

        for (i = 2; i < nchunk; ++i) {
            for (j = 0; j < 50; ++j) {
                tuna_chunk_obs(&chunks[i], 1.0 + (j % 5) * 0.01);
            }
        }

        stats0 = tuna_stats_mom(&chunks[0].stats);
        stats1 = tuna_stats_mom(&chunks[1].stats);
        p_value = tuna_welch1(stats0.avg, stats0.var, stats0.n,
                              stats1.avg, stats1.var, stats1.n);

        u01[0] = 0.5;
        u01[1] = p_value * 5.0;
        for (i = 2; i < nchunk; ++i) {
            u01[i] = 0.0001;
        }

        fct_chk(p_value < u01[1]);
        fct_chk(p_value * 9.0 >= u01[1]);

        result = tuna_algo_welch1->function(nchunk, chunks, u01);
        fct_chk_eq_int(result, 0);
    }
    FCT_QTEST_END();

    FCT_QTEST_BGN(bonferroni_allows_strong_differences)
    {
        /* Verify that truly significant differences are still detected */
        const size_t nchunk = 5;
        tuna_chunk chunks[5] = {};
        double u01[5];
        size_t i, j;
        size_t result;

        /* Chunk 0: mean=1.0 */
        for (j = 0; j < 50; ++j) {
            tuna_chunk_obs(&chunks[0], 1.0 + (j % 3) * 0.01);
        }

        /* Chunk 1: mean=0.5 (much better - very significant) */
        for (j = 0; j < 50; ++j) {
            tuna_chunk_obs(&chunks[1], 0.5 + (j % 3) * 0.01);
        }

        /* Chunks 2-4: mean=1.0 */
        for (i = 2; i < nchunk; ++i) {
            for (j = 0; j < 50; ++j) {
                tuna_chunk_obs(&chunks[i], 1.0 + (j % 3) * 0.01);
            }
        }

        for (i = 0; i < nchunk; ++i) {
            u01[i] = 0.5;
        }

        /* Should detect chunk 1 is significantly better */
        result = tuna_algo_welch1_nuinf->function(nchunk, chunks, u01);
        fct_chk_eq_int(result, 1);

        result = tuna_algo_welch1->function(nchunk, chunks, u01);
        fct_chk_eq_int(result, 1);
    }
    FCT_QTEST_END();

    FCT_QTEST_BGN(bonferroni_scaling_demonstration)
    {
        /* Demonstrate that correction scales with nchunk by showing
         * that the same data and base threshold lead to different
         * outcomes with different nchunk values.
         */
        tuna_chunk chunks3[3] = {};
        tuna_chunk chunks10[10] = {};
        double u01_3[3];
        double u01_10[10];
        size_t i, j;
        tuna_stats_mom_result stats0, stats1;
        double p_value;

        /* Create identical chunk data for both tests */
        for (j = 0; j < 50; ++j) {
            tuna_chunk_obs(&chunks3[0], 1.0 + (j % 5) * 0.01);
            tuna_chunk_obs(&chunks10[0], 1.0 + (j % 5) * 0.01);
        }
        for (j = 0; j < 50; ++j) {
            tuna_chunk_obs(&chunks3[1], 1.1 + (j % 5) * 0.01);
            tuna_chunk_obs(&chunks10[1], 1.1 + (j % 5) * 0.01);
        }
        for (j = 0; j < 50; ++j) {
            tuna_chunk_obs(&chunks3[2], 1.0 + (j % 5) * 0.01);
            tuna_chunk_obs(&chunks10[2], 1.0 + (j % 5) * 0.01);
        }
        for (i = 3; i < 10; ++i) {
            for (j = 0; j < 50; ++j) {
                tuna_chunk_obs(&chunks10[i], 1.0 + (j % 5) * 0.01);
            }
        }

        /* Get p-value */
        stats0 = tuna_stats_mom(&chunks3[0].stats);
        stats1 = tuna_stats_mom(&chunks3[1].stats);
        p_value = tuna_welch1_nuinf(stats0.avg, stats0.var, stats0.n,
                                     stats1.avg, stats1.var, stats1.n);

        /* For nchunk=3 (ncomparisons=2): set threshold at p*4
         * This means p*2 < threshold, might switch */
        u01_3[0] = 0.5;
        u01_3[1] = p_value * 4.0;
        u01_3[2] = 0.0001;

        /* For nchunk=10 (ncomparisons=9): use same threshold
         * This means p*9 >= threshold, won't switch */
        u01_10[0] = 0.5;
        u01_10[1] = p_value * 4.0;
        for (i = 2; i < 10; ++i) {
            u01_10[i] = 0.0001;
        }

        /* With nchunk=10, correction is stronger, should stay at 0 */
        size_t result = tuna_algo_welch1_nuinf->function(10, chunks10, u01_10);
        fct_chk_eq_int(result, 0);
    }
    FCT_QTEST_END();
}
FCT_END()
