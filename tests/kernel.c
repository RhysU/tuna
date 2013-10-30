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
    }
    FCT_QTEST_END();

}
FCT_END()
