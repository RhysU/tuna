//--------------------------------------------------------------------------
//
// Copyright (C) 2013 Rhys Ulerich
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
    FCT_QTEST_BGN(extern_declarations)
    {
        // Test that the extern algorithm handles are accessible
        fct_chk(tuna_algo_welch1 != NULL);
        fct_chk(tuna_algo_welch1_nuinf != NULL);
        fct_chk(tuna_algo_zero != NULL);

        // Test that they are distinct
        fct_chk(tuna_algo_welch1 != tuna_algo_welch1_nuinf);
        fct_chk(tuna_algo_welch1 != tuna_algo_zero);
        fct_chk(tuna_algo_welch1_nuinf != tuna_algo_zero);
    }
    FCT_QTEST_END();

    FCT_QTEST_BGN(algorithm_names)
    {
        // Test retrieving names from algorithm handles
        fct_chk_eq_str(tuna_algo_name(tuna_algo_welch1), "welch1");
        fct_chk_eq_str(tuna_algo_name(tuna_algo_welch1_nuinf), "welch1_nuinf");
        fct_chk_eq_str(tuna_algo_name(tuna_algo_zero), "zero");

        // Test NULL algorithm returns "unknown"
        fct_chk_eq_str(tuna_algo_name(NULL), "unknown");

        // Test accessing name directly from struct
        fct_chk_eq_str(tuna_algo_welch1->name, "welch1");
        fct_chk_eq_str(tuna_algo_welch1_nuinf->name, "welch1_nuinf");
        fct_chk_eq_str(tuna_algo_zero->name, "zero");
    }
    FCT_QTEST_END();

    FCT_QTEST_BGN(invoke_zero_algorithm)
    {
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

    FCT_QTEST_BGN(invoke_welch1_algorithm)
    {
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

    FCT_QTEST_BGN(invoke_welch1_nuinf_algorithm)
    {
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

    FCT_QTEST_BGN(tuna_algo_default_behavior)
    {
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

    FCT_QTEST_BGN(algorithm_function_pointer_validity)
    {
        // Test that function pointers are valid
        fct_chk(tuna_algo_welch1->function != NULL);
        fct_chk(tuna_algo_welch1_nuinf->function != NULL);
        fct_chk(tuna_algo_zero->function != NULL);

        // Test that they are distinct functions
        fct_chk(tuna_algo_welch1->function != tuna_algo_zero->function);
        fct_chk(tuna_algo_welch1_nuinf->function != tuna_algo_zero->function);
        fct_chk(tuna_algo_welch1->function != tuna_algo_welch1_nuinf->function);
    }
    FCT_QTEST_END();

    FCT_QTEST_BGN(algorithm_in_tuna_site)
    {
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
    }
    FCT_QTEST_END();

}
FCT_END()
