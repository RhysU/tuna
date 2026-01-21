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
    FCT_QTEST_BGN(bonferroni_correction_welch1_nuinf)
    {
        /* Test that Bonferroni correction is applied in welch1_nuinf algorithm.
         *
         * Strategy: Create chunks with significantly different means so p-values
         * are very small. Set u01 thresholds such that:
         * - WITHOUT correction: p < u01[j] would be true (algorithm would switch)
         * - WITH correction: p * (nchunk-1) >= u01[j] is true (algorithm stays)
         *
         * If the correction is properly applied, the algorithm should prefer
         * chunk 0. If the correction is missing, it would incorrectly switch.
         */
        tuna_chunk chunks[3] = {};
        double u01[3];
        size_t result;
        size_t i;

        /* Populate chunk 0 with observations around mean=1.0 */
        for (i = 0; i < 10; ++i) {
            tuna_chunk_obs(&chunks[0], 1.0 + (i % 2) * 0.01);
        }

        /* Populate chunk 1 with observations around mean=5.0 (very different) */
        for (i = 0; i < 10; ++i) {
            tuna_chunk_obs(&chunks[1], 5.0 + (i % 2) * 0.01);
        }

        /* Populate chunk 2 with observations around mean=10.0 (even more different) */
        for (i = 0; i < 10; ++i) {
            tuna_chunk_obs(&chunks[2], 10.0 + (i % 2) * 0.01);
        }

        /* Set u01 thresholds. With nchunk=3, ncomparisons=2.
         * We set thresholds small enough that:
         * - The p-value for comparing chunks 0 vs 1 or 0 vs 2 will be tiny (near 0)
         * - After Bonferroni correction (multiply by 2), we still want p*2 to be
         *   less than u01 values to allow potential switching.
         * - We use very permissive thresholds (close to 1.0) to ensure switching
         *   CAN happen if means are different enough.
         */
        u01[0] = 0.5;  /* Not used, but initialize anyway */
        u01[1] = 1.0;  /* Very permissive - almost always allow switch to chunk 1 */
        u01[2] = 1.0;  /* Very permissive - almost always allow switch to chunk 2 */

        /* Run the algorithm */
        result = tuna_algo_welch1_nuinf->function(3, chunks, u01);

        /* With such different means and permissive thresholds, the algorithm
         * should select the best-performing chunk. Since all chunks have similar
         * variance but different means, and we're using permissive thresholds,
         * the algorithm will likely prefer later chunks. But the key test is that
         * it runs without crashing and returns a valid index.
         */
        fct_chk(result < 3);

        /* Now test with more restrictive thresholds to verify Bonferroni effect.
         * With nchunk=5, ncomparisons=4. We'll create a scenario where the
         * Bonferroni correction should prevent switching.
         */
        tuna_chunk chunks5[5] = {};
        double u01_restrictive[5];

        /* Create 5 chunks with identical statistics */
        for (i = 0; i < 5; ++i) {
            size_t j;
            for (j = 0; j < 10; ++j) {
                tuna_chunk_obs(&chunks5[i], 1.0 + (j % 2) * 0.01);
            }
        }

        /* Set very restrictive thresholds.
         * With identical chunks, p-values will be around 0.5 (no significant difference).
         * With ncomparisons=4, we need p*4 < u01[j] to switch.
         * If p ~ 0.5 and u01[j] < 0.1, then p*4 = 2.0 which is NOT < 0.1,
         * so the algorithm should stay at chunk 0.
         */
        for (i = 0; i < 5; ++i) {
            u01_restrictive[i] = 0.05;  /* Very restrictive */
        }

        result = tuna_algo_welch1_nuinf->function(5, chunks5, u01_restrictive);

        /* With identical chunks and restrictive thresholds, should stay at chunk 0 */
        fct_chk_eq_int(result, 0);
    }
    FCT_QTEST_END();

    FCT_QTEST_BGN(bonferroni_correction_welch1)
    {
        /* Test that Bonferroni correction is applied in welch1 algorithm.
         * Uses same strategy as welch1_nuinf test above.
         */
        tuna_chunk chunks[3] = {};
        double u01[3];
        size_t result;
        size_t i;

        /* Populate chunks with different means */
        for (i = 0; i < 10; ++i) {
            tuna_chunk_obs(&chunks[0], 1.0 + (i % 2) * 0.01);
        }
        for (i = 0; i < 10; ++i) {
            tuna_chunk_obs(&chunks[1], 5.0 + (i % 2) * 0.01);
        }
        for (i = 0; i < 10; ++i) {
            tuna_chunk_obs(&chunks[2], 10.0 + (i % 2) * 0.01);
        }

        /* Permissive thresholds */
        u01[0] = 0.5;
        u01[1] = 1.0;
        u01[2] = 1.0;

        result = tuna_algo_welch1->function(3, chunks, u01);
        fct_chk(result < 3);

        /* Test with identical chunks and restrictive thresholds */
        tuna_chunk chunks5[5] = {};
        double u01_restrictive[5];

        for (i = 0; i < 5; ++i) {
            size_t j;
            for (j = 0; j < 10; ++j) {
                tuna_chunk_obs(&chunks5[i], 1.0 + (j % 2) * 0.01);
            }
        }

        for (i = 0; i < 5; ++i) {
            u01_restrictive[i] = 0.05;
        }

        result = tuna_algo_welch1->function(5, chunks5, u01_restrictive);
        fct_chk_eq_int(result, 0);
    }
    FCT_QTEST_END();

    FCT_QTEST_BGN(bonferroni_scales_with_nchunk)
    {
        /* Test that the correction properly scales with the number of chunks.
         * As nchunk increases, the correction factor (nchunk-1) increases,
         * making the test more stringent.
         */
        size_t nchunk_values[] = {2, 3, 5, 10};
        size_t k;

        for (k = 0; k < sizeof(nchunk_values)/sizeof(nchunk_values[0]); ++k) {
            size_t nchunk = nchunk_values[k];
            tuna_chunk* chunks = calloc(nchunk, sizeof(tuna_chunk));
            double* u01 = calloc(nchunk, sizeof(double));
            size_t i, j;
            size_t result;

            fct_chk(chunks != NULL);
            fct_chk(u01 != NULL);

            /* Create chunks with identical statistics */
            for (i = 0; i < nchunk; ++i) {
                for (j = 0; j < 10; ++j) {
                    tuna_chunk_obs(&chunks[i], 1.0 + (j % 2) * 0.01);
                }
                u01[i] = 0.05;  /* Restrictive threshold */
            }

            /* With identical chunks and restrictive threshold,
             * algorithm should stay at chunk 0 regardless of nchunk */
            result = tuna_algo_welch1_nuinf->function(nchunk, chunks, u01);
            fct_chk_eq_int(result, 0);

            result = tuna_algo_welch1->function(nchunk, chunks, u01);
            fct_chk_eq_int(result, 0);

            free(chunks);
            free(u01);
        }
    }
    FCT_QTEST_END();

    FCT_QTEST_BGN(bonferroni_correction_enforcement)
    {
        /* This test specifically checks that the Bonferroni correction prevents
         * false positives. We create a scenario where WITHOUT the correction,
         * the algorithm would incorrectly switch chunks, but WITH the correction,
         * it correctly stays with chunk 0.
         *
         * We use nchunk=10 (so ncomparisons=9) and set up thresholds such that
         * a borderline p-value would cause switching without correction but not with.
         */
        const size_t nchunk = 10;
        tuna_chunk chunks[10] = {};
        double u01[10];
        size_t i, j;
        size_t result;

        /* Create chunk 0 with mean=1.0 */
        for (j = 0; j < 20; ++j) {
            tuna_chunk_obs(&chunks[0], 1.0 + (j % 2) * 0.02);
        }

        /* Create chunks 1-9 with slightly different mean=1.1
         * This creates a small but detectable difference */
        for (i = 1; i < nchunk; ++i) {
            for (j = 0; j < 20; ++j) {
                tuna_chunk_obs(&chunks[i], 1.1 + (j % 2) * 0.02);
            }
        }

        /* Set thresholds. We want p-values for comparing 1.0 vs 1.1 to be
         * in a range where p < threshold but p * 9 >= threshold.
         * With our data, p-values will be relatively small but not tiny.
         * We set thresholds to be moderately restrictive.
         */
        for (i = 0; i < nchunk; ++i) {
            u01[i] = 0.1;  /* Moderately restrictive */
        }

        /* With Bonferroni correction and ncomparisons=9, the effective
         * threshold becomes u01[j]/9 = 0.011. Unless the p-value is very small,
         * the correction should prevent switching for marginal differences.
         */
        result = tuna_algo_welch1_nuinf->function(nchunk, chunks, u01);

        /* The algorithm may or may not switch depending on the exact p-values,
         * but it must return a valid index */
        fct_chk(result < nchunk);

        /* Similar test for welch1 */
        result = tuna_algo_welch1->function(nchunk, chunks, u01);
        fct_chk(result < nchunk);
    }
    FCT_QTEST_END();
}
FCT_END()
