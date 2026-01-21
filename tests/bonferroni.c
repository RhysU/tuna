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
    // This test is designed to FAIL if Bonferroni correction is removed.
    //
    // Strategy: Create chunks where chunk 0 is better (lower mean) but chunk 1
    // is worse (higher mean). Calculate the actual p-value and set u01[1] such
    // that p < u01[1] BUT p*ncomparisons >= u01[1]. This creates a boundary
    // condition that distinguishes between corrected and uncorrected behavior.
    FCT_QTEST_BGN(bonferroni_without_correction_would_switch)
    {
        const size_t nchunk = 10;  /* ncomparisons = 9 */
        tuna_chunk chunks[10] = {};
        double u01[10];
        tuna_stats_mom_result stats0, stats1;
        double p_value;
        size_t i, j;

        // Populate chunks: 0=better (1.0), 1=worse (1.1), rest same as 0
        for (j = 0; j < 50; ++j) {
            double val0 = 1.0 + (j % 5) * 0.01;
            double val1 = 1.1 + (j % 5) * 0.01;
            tuna_chunk_obs(&chunks[0], val0);
            tuna_chunk_obs(&chunks[1], val1);
            for (i = 2; i < nchunk; ++i) {
                tuna_chunk_obs(&chunks[i], val0);
            }
        }

        // Calculate p-value and design threshold to catch Bonferroni effect
        stats0 = tuna_stats_mom(&chunks[0].stats);
        stats1 = tuna_stats_mom(&chunks[1].stats);
        p_value = tuna_welch1_nuinf(stats0.avg, stats0.var, stats0.n,
                                     stats1.avg, stats1.var, stats1.n);

        u01[0] = 0.5;
        u01[1] = p_value * 5.0;  // Between p and p*9
        for (i = 2; i < nchunk; ++i) {
            u01[i] = 0.0001;
        }

        // Verify threshold design
        fct_chk(p_value < u01[1]);         // Would switch without correction
        fct_chk(p_value * 9.0 >= u01[1]);  // Stays with correction

        // With Bonferroni correction, should stay at chunk 0
        size_t result = tuna_algo_welch1_nuinf->function(nchunk, chunks, u01);
        fct_chk_eq_int(result, 0);
    }
    FCT_QTEST_END();

    // Same test for standard welch1 algorithm
    FCT_QTEST_BGN(bonferroni_welch1_variant)
    {
        const size_t nchunk = 10;
        tuna_chunk chunks[10] = {};
        double u01[10];
        tuna_stats_mom_result stats0, stats1;
        double p_value;
        size_t i, j;

        for (j = 0; j < 50; ++j) {
            double val0 = 1.0 + (j % 5) * 0.01;
            double val1 = 1.1 + (j % 5) * 0.01;
            tuna_chunk_obs(&chunks[0], val0);
            tuna_chunk_obs(&chunks[1], val1);
            for (i = 2; i < nchunk; ++i) {
                tuna_chunk_obs(&chunks[i], val0);
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

        size_t result = tuna_algo_welch1->function(nchunk, chunks, u01);
        fct_chk_eq_int(result, 0);
    }
    FCT_QTEST_END();

    // Verify truly significant differences are still detected
    FCT_QTEST_BGN(bonferroni_allows_strong_differences)
    {
        const size_t nchunk = 5;
        tuna_chunk chunks[5] = {};
        double u01[5] = {0.5, 0.5, 0.5, 0.5, 0.5};
        size_t i, j;

        for (j = 0; j < 50; ++j) {
            double val_normal = 1.0 + (j % 3) * 0.01;
            double val_better = 0.5 + (j % 3) * 0.01;
            tuna_chunk_obs(&chunks[0], val_normal);
            tuna_chunk_obs(&chunks[1], val_better);  // Much better
            for (i = 2; i < nchunk; ++i) {
                tuna_chunk_obs(&chunks[i], val_normal);
            }
        }

        // Should detect chunk 1 is significantly better
        size_t result = tuna_algo_welch1_nuinf->function(nchunk, chunks, u01);
        fct_chk_eq_int(result, 1);

        result = tuna_algo_welch1->function(nchunk, chunks, u01);
        fct_chk_eq_int(result, 1);
    }
    FCT_QTEST_END();

    FCT_QTEST_BGN(bonferroni_scaling_demonstration)
    {
        // Demonstrate correction scales with nchunk
        tuna_chunk chunks10[10] = {};
        double u01_10[10];
        tuna_stats_mom_result stats0, stats1;
        double p_value;
        size_t i, j;

        for (j = 0; j < 50; ++j) {
            double val0 = 1.0 + (j % 5) * 0.01;
            double val1 = 1.1 + (j % 5) * 0.01;
            tuna_chunk_obs(&chunks10[0], val0);
            tuna_chunk_obs(&chunks10[1], val1);
            for (i = 2; i < 10; ++i) {
                tuna_chunk_obs(&chunks10[i], val0);
            }
        }

        stats0 = tuna_stats_mom(&chunks10[0].stats);
        stats1 = tuna_stats_mom(&chunks10[1].stats);
        p_value = tuna_welch1_nuinf(stats0.avg, stats0.var, stats0.n,
                                     stats1.avg, stats1.var, stats1.n);

        // For nchunk=10 (ncomparisons=9)
        // threshold at p*4 means p*9 >= threshold
        u01_10[0] = 0.5;
        u01_10[1] = p_value * 4.0;
        for (i = 2; i < 10; ++i) {
            u01_10[i] = 0.0001;
        }

        // With nchunk=10, correction is stronger, should stay at 0
        size_t result = tuna_algo_welch1_nuinf->function(10, chunks10, u01_10);
        fct_chk_eq_int(result, 0);
    }
    FCT_QTEST_END();
}
FCT_END()
