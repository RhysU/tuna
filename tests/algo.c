//--------------------------------------------------------------------------
//
// Copyright (C) 2013, 2026 Rhys Ulerich
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
    FCT_QTEST_BGN(extern_declarations) {
        // Test that the extern algorithm handles are accessible
        fct_chk(tuna_algo_welch1 != NULL);
        fct_chk(tuna_algo_welch1_nuinf != NULL);
        fct_chk(tuna_algo_zero != NULL);
        fct_chk(tuna_algo_uniform != NULL);
        fct_chk(tuna_algo_thompson != NULL);
        fct_chk(tuna_algo_ucb1 != NULL);

        // Test that they are distinct
        fct_chk(tuna_algo_welch1 != tuna_algo_welch1_nuinf);
        fct_chk(tuna_algo_welch1 != tuna_algo_zero);
        fct_chk(tuna_algo_welch1 != tuna_algo_uniform);
        fct_chk(tuna_algo_welch1 != tuna_algo_thompson);
        fct_chk(tuna_algo_welch1 != tuna_algo_ucb1);
        fct_chk(tuna_algo_welch1_nuinf != tuna_algo_zero);
        fct_chk(tuna_algo_welch1_nuinf != tuna_algo_uniform);
        fct_chk(tuna_algo_welch1_nuinf != tuna_algo_thompson);
        fct_chk(tuna_algo_welch1_nuinf != tuna_algo_ucb1);
        fct_chk(tuna_algo_zero != tuna_algo_uniform);
        fct_chk(tuna_algo_zero != tuna_algo_thompson);
        fct_chk(tuna_algo_zero != tuna_algo_ucb1);
        fct_chk(tuna_algo_uniform != tuna_algo_thompson);
        fct_chk(tuna_algo_uniform != tuna_algo_ucb1);
        fct_chk(tuna_algo_thompson != tuna_algo_ucb1);
    }
    FCT_QTEST_END();

    FCT_QTEST_BGN(algorithm_names) {
        // Test retrieving names from algorithm handles
        fct_chk_eq_str(tuna_algo_name(tuna_algo_welch1), "welch1");
        fct_chk_eq_str(tuna_algo_name(tuna_algo_welch1_nuinf), "welch1_nuinf");
        fct_chk_eq_str(tuna_algo_name(tuna_algo_zero), "zero");
        fct_chk_eq_str(tuna_algo_name(tuna_algo_uniform), "uniform");
        fct_chk_eq_str(tuna_algo_name(tuna_algo_thompson), "thompson");
        fct_chk_eq_str(tuna_algo_name(tuna_algo_ucb1), "ucb1");

        // Test NULL algorithm returns "unknown"
        fct_chk_eq_str(tuna_algo_name(NULL), "unknown");

        // Test accessing name directly from struct
        fct_chk_eq_str(tuna_algo_welch1->name, "welch1");
        fct_chk_eq_str(tuna_algo_welch1_nuinf->name, "welch1_nuinf");
        fct_chk_eq_str(tuna_algo_zero->name, "zero");
        fct_chk_eq_str(tuna_algo_uniform->name, "uniform");
        fct_chk_eq_str(tuna_algo_thompson->name, "thompson");
        fct_chk_eq_str(tuna_algo_ucb1->name, "ucb1");
    }
    FCT_QTEST_END();

    FCT_QTEST_BGN(invoke_zero_algorithm) {
        // Test invoking the zero algorithm function
        tuna_chunk chunks[3] = {};
        double u01[3] = {0.5, 0.5, 0.5};

        // zero algorithm should always return index 0
        size_t result = tuna_algo_zero->function(3, chunks, u01);
        fct_chk_eq_int(result, 0);

        // Test with different inputs
        result = tuna_algo_zero->function(5, chunks, u01);
        fct_chk_eq_int(result, 0);
    }
    FCT_QTEST_END();

    FCT_QTEST_BGN(invoke_welch1_algorithm) {
        // Test invoking welch1 algorithm with mock data
        tuna_chunk chunks[2] = {};
        double u01[2] = {0.5, 0.5};

        // Add observations to first chunk
        tuna_chunk_obs(&chunks[0], 1.0);
        tuna_chunk_obs(&chunks[0], 1.0);
        tuna_chunk_obs(&chunks[0], 1.0);
        tuna_chunk_obs(&chunks[0], 1.0);

        // Second chunk starts with insufficient data
        size_t result = tuna_algo_welch1->function(2, chunks, u01);
        // Should return 0 or 1 depending on statistics
        fct_chk(result < 2);

        // Add observations to second chunk
        tuna_chunk_obs(&chunks[1], 2.0);
        tuna_chunk_obs(&chunks[1], 2.0);
        tuna_chunk_obs(&chunks[1], 2.0);
        tuna_chunk_obs(&chunks[1], 2.0);

        result = tuna_algo_welch1->function(2, chunks, u01);
        fct_chk(result < 2);
    }
    FCT_QTEST_END();

    FCT_QTEST_BGN(invoke_welch1_nuinf_algorithm) {
        // Test invoking welch1_nuinf algorithm with mock data
        tuna_chunk chunks[2] = {};
        double u01[2] = {0.5, 0.5};

        // Add observations to chunks
        tuna_chunk_obs(&chunks[0], 1.0);
        tuna_chunk_obs(&chunks[0], 1.0);
        tuna_chunk_obs(&chunks[0], 1.0);
        tuna_chunk_obs(&chunks[0], 1.0);

        tuna_chunk_obs(&chunks[1], 2.0);
        tuna_chunk_obs(&chunks[1], 2.0);
        tuna_chunk_obs(&chunks[1], 2.0);
        tuna_chunk_obs(&chunks[1], 2.0);

        size_t result = tuna_algo_welch1_nuinf->function(2, chunks, u01);
        fct_chk(result < 2);
    }
    FCT_QTEST_END();

    FCT_QTEST_BGN(tuna_algo_default_behavior) {
        // Test default algorithm selection
        const tuna_algo* algo;

        // When nchunk < 2, should return zero algorithm
        algo = tuna_algo_default(0);
        fct_chk(algo == tuna_algo_zero);

        algo = tuna_algo_default(1);
        fct_chk(algo == tuna_algo_zero);

        // When nchunk >= 2, should return welch1 (default)
        algo = tuna_algo_default(2);
        fct_chk(algo == tuna_algo_welch1);

        algo = tuna_algo_default(10);
        fct_chk(algo == tuna_algo_welch1);
    }
    FCT_QTEST_END();

    FCT_QTEST_BGN(algorithm_function_pointer_validity) {
        // Test that function pointers are valid
        fct_chk(tuna_algo_welch1->function != NULL);
        fct_chk(tuna_algo_welch1_nuinf->function != NULL);
        fct_chk(tuna_algo_zero->function != NULL);
        fct_chk(tuna_algo_uniform->function != NULL);
        fct_chk(tuna_algo_thompson->function != NULL);
        fct_chk(tuna_algo_ucb1->function != NULL);

        // Test that they are distinct functions
        fct_chk(tuna_algo_welch1->function != tuna_algo_zero->function);
        fct_chk(tuna_algo_welch1->function != tuna_algo_uniform->function);
        fct_chk(tuna_algo_welch1->function != tuna_algo_thompson->function);
        fct_chk(tuna_algo_welch1->function != tuna_algo_ucb1->function);
        fct_chk(tuna_algo_welch1_nuinf->function != tuna_algo_zero->function);
        fct_chk(tuna_algo_welch1_nuinf->function != tuna_algo_uniform->function);
        fct_chk(tuna_algo_welch1_nuinf->function != tuna_algo_thompson->function);
        fct_chk(tuna_algo_welch1_nuinf->function != tuna_algo_ucb1->function);
        fct_chk(tuna_algo_welch1->function != tuna_algo_welch1_nuinf->function);
        fct_chk(tuna_algo_zero->function != tuna_algo_uniform->function);
        fct_chk(tuna_algo_zero->function != tuna_algo_thompson->function);
        fct_chk(tuna_algo_zero->function != tuna_algo_ucb1->function);
        fct_chk(tuna_algo_uniform->function != tuna_algo_thompson->function);
        fct_chk(tuna_algo_uniform->function != tuna_algo_ucb1->function);
        fct_chk(tuna_algo_thompson->function != tuna_algo_ucb1->function);
    }
    FCT_QTEST_END();

    FCT_QTEST_BGN(algorithm_in_tuna_site) {
        // Test using algorithm handles with tuna_site
        tuna_site site = {};

        site.algo = tuna_algo_zero;
        fct_chk(site.algo == tuna_algo_zero);
        fct_chk_eq_str(tuna_algo_name(site.algo), "zero");

        site.algo = tuna_algo_welch1;
        fct_chk(site.algo == tuna_algo_welch1);
        fct_chk_eq_str(tuna_algo_name(site.algo), "welch1");

        site.algo = tuna_algo_welch1_nuinf;
        fct_chk(site.algo == tuna_algo_welch1_nuinf);
        fct_chk_eq_str(tuna_algo_name(site.algo), "welch1_nuinf");

        site.algo = tuna_algo_uniform;
        fct_chk(site.algo == tuna_algo_uniform);
        fct_chk_eq_str(tuna_algo_name(site.algo), "uniform");

        site.algo = tuna_algo_thompson;
        fct_chk(site.algo == tuna_algo_thompson);
        fct_chk_eq_str(tuna_algo_name(site.algo), "thompson");

        site.algo = tuna_algo_ucb1;
        fct_chk(site.algo == tuna_algo_ucb1);
        fct_chk_eq_str(tuna_algo_name(site.algo), "ucb1");
    }
    FCT_QTEST_END();

    FCT_QTEST_BGN(invoke_uniform_algorithm) {
        // Test invoking uniform algorithm
        tuna_chunk chunks[5] = {};
        double u01[5];
        size_t result;
        int i;

        // Test uniform selection with various random values
        // u01[0] = 0.0 should select chunk 0
        u01[0] = 0.0;
        for (i = 1; i < 5; ++i) {
            u01[i] = 0.5;
        }
        result = tuna_algo_uniform->function(5, chunks, u01);
        fct_chk_eq_int(result, 0);

        // u01[0] = 0.199 should select chunk 0 (0.199 * 5 = 0.995 -> floor = 0)
        u01[0] = 0.199;
        result = tuna_algo_uniform->function(5, chunks, u01);
        fct_chk_eq_int(result, 0);

        // u01[0] = 0.2 should select chunk 1 (0.2 * 5 = 1.0)
        u01[0] = 0.2;
        result = tuna_algo_uniform->function(5, chunks, u01);
        fct_chk_eq_int(result, 1);

        // u01[0] = 0.4 should select chunk 2 (0.4 * 5 = 2.0)
        u01[0] = 0.4;
        result = tuna_algo_uniform->function(5, chunks, u01);
        fct_chk_eq_int(result, 2);

        // u01[0] = 0.6 should select chunk 3 (0.6 * 5 = 3.0)
        u01[0] = 0.6;
        result = tuna_algo_uniform->function(5, chunks, u01);
        fct_chk_eq_int(result, 3);

        // u01[0] = 0.8 should select chunk 4 (0.8 * 5 = 4.0)
        u01[0] = 0.8;
        result = tuna_algo_uniform->function(5, chunks, u01);
        fct_chk_eq_int(result, 4);

        // u01[0] close to 1.0 should select chunk 4 (0.999 * 5 = 4.995 -> floor = 4)
        u01[0] = 0.999;
        result = tuna_algo_uniform->function(5, chunks, u01);
        fct_chk_eq_int(result, 4);

        // Test with different number of chunks
        u01[0] = 0.5;
        result = tuna_algo_uniform->function(2, chunks, u01);
        fct_chk_eq_int(result, 1);

        u01[0] = 0.49;
        result = tuna_algo_uniform->function(2, chunks, u01);
        fct_chk_eq_int(result, 0);

        // Test that uniform ignores chunk statistics (unlike welch algorithms)
        tuna_chunk_obs(&chunks[0], 100.0);
        tuna_chunk_obs(&chunks[0], 100.0);
        tuna_chunk_obs(&chunks[1], 1.0);
        tuna_chunk_obs(&chunks[1], 1.0);

        u01[0] = 0.1;
        result = tuna_algo_uniform->function(2, chunks, u01);
        fct_chk_eq_int(result, 0);  // Should still pick 0 based on u01[0]
    }
    FCT_QTEST_END();

    FCT_QTEST_BGN(invoke_thompson_algorithm) {
        // Test invoking Thompson Sampling algorithm
        // Note: tuna_chunk uses 3 burn-in observations before stats accumulate,
        // so we need 5+ observations to get stats.n >= 2
        tuna_chunk chunks[3] = {};
        double u01[3];
        size_t result;
        int i;

        // With no data, should return first chunk needing data
        u01[0] = 0.5; u01[1] = 0.5; u01[2] = 0.5;
        result = tuna_algo_thompson->function(3, chunks, u01);
        fct_chk_eq_int(result, 0);

        // Add enough observations to first chunk (high cost) for stats.n >= 2
        // 3 burn-in + 2 actual = 5 observations needed
        for (i = 0; i < 6; ++i) {
            tuna_chunk_obs(&chunks[0], 10.0);
        }

        // Second chunk still needs data
        result = tuna_algo_thompson->function(3, chunks, u01);
        fct_chk_eq_int(result, 1);

        // Add observations to second chunk (low cost)
        for (i = 0; i < 6; ++i) {
            tuna_chunk_obs(&chunks[1], 1.0);
        }

        // Third chunk still needs data
        result = tuna_algo_thompson->function(3, chunks, u01);
        fct_chk_eq_int(result, 2);

        // Add observations to third chunk (medium cost)
        for (i = 0; i < 6; ++i) {
            tuna_chunk_obs(&chunks[2], 5.0);
        }

        // Now with all data, should return valid index
        result = tuna_algo_thompson->function(3, chunks, u01);
        fct_chk(result < 3);

        // With u01 values that give extreme z-scores,
        // low u01 -> negative z -> should favor low mean chunks
        u01[0] = 0.01; u01[1] = 0.01; u01[2] = 0.01;
        result = tuna_algo_thompson->function(3, chunks, u01);
        fct_chk(result < 3);

        // High u01 -> positive z -> increases sampled values
        u01[0] = 0.99; u01[1] = 0.99; u01[2] = 0.99;
        result = tuna_algo_thompson->function(3, chunks, u01);
        fct_chk(result < 3);
    }
    FCT_QTEST_END();

    FCT_QTEST_BGN(invoke_ucb1_algorithm) {
        // Test invoking UCB1 algorithm
        // Note: tuna_chunk uses 3 burn-in observations before stats accumulate,
        // so we need 5+ observations to get stats.n >= 2
        tuna_chunk chunks[3] = {};
        double u01[3] = {0.5, 0.5, 0.5};
        size_t result, result2;
        int i;

        // With no data, should return first chunk
        result = tuna_algo_ucb1->function(3, chunks, u01);
        fct_chk_eq_int(result, 0);

        // Add enough observations to first chunk (high cost) for stats.n >= 2
        for (i = 0; i < 6; ++i) {
            tuna_chunk_obs(&chunks[0], 10.0);
        }

        // Second chunk needs data
        result = tuna_algo_ucb1->function(3, chunks, u01);
        fct_chk_eq_int(result, 1);

        // Add observations to second chunk (low cost)
        for (i = 0; i < 6; ++i) {
            tuna_chunk_obs(&chunks[1], 1.0);
        }

        // Third chunk needs data
        result = tuna_algo_ucb1->function(3, chunks, u01);
        fct_chk_eq_int(result, 2);

        // Add observations to third chunk (medium cost)
        for (i = 0; i < 6; ++i) {
            tuna_chunk_obs(&chunks[2], 5.0);
        }

        // Now with all data, UCB1 should return valid index
        // UCB1 is deterministic (ignores u01), so result should be consistent
        result = tuna_algo_ucb1->function(3, chunks, u01);
        fct_chk(result < 3);

        // Run again with different u01 - should get same result (deterministic)
        u01[0] = 0.1; u01[1] = 0.9; u01[2] = 0.5;
        result2 = tuna_algo_ucb1->function(3, chunks, u01);
        fct_chk_eq_int(result, result2);

        // With significant difference in means and equal samples,
        // UCB1 should favor the lower cost chunk
        // chunk[1] has mean=1.0, chunk[0] has mean=10.0, chunk[2] has mean=5.0
        // LCB = mean - exploration_bonus, so chunk[1] should have lowest LCB
        fct_chk_eq_int(result, 1);
    }
    FCT_QTEST_END();

    FCT_QTEST_BGN(thompson_vs_ucb1_determinism) {
        // Test that UCB1 is deterministic while Thompson varies with u01
        // Note: tuna_chunk uses 3 burn-in observations before stats accumulate
        tuna_chunk chunks[2] = {};
        double u01_low[2] = {0.1, 0.1};
        double u01_high[2] = {0.9, 0.9};
        size_t result_ucb1_low, result_ucb1_high;
        size_t result_thompson_low, result_thompson_high;
        int i;

        // Add enough observations to both chunks (6 = 3 burn-in + 3 actual)
        for (i = 0; i < 6; ++i) {
            tuna_chunk_obs(&chunks[0], 2.0);
        }

        for (i = 0; i < 6; ++i) {
            tuna_chunk_obs(&chunks[1], 2.1);
        }

        // UCB1 should give same result regardless of u01
        result_ucb1_low = tuna_algo_ucb1->function(2, chunks, u01_low);
        result_ucb1_high = tuna_algo_ucb1->function(2, chunks, u01_high);
        fct_chk_eq_int(result_ucb1_low, result_ucb1_high);

        // Thompson may give different results based on u01
        // (though not guaranteed for every specific u01 pair)
        result_thompson_low = tuna_algo_thompson->function(2, chunks, u01_low);
        result_thompson_high = tuna_algo_thompson->function(2, chunks, u01_high);
        fct_chk(result_thompson_low < 2);
        fct_chk(result_thompson_high < 2);
    }
    FCT_QTEST_END();

}
FCT_END()
