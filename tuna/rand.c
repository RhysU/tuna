/*
 * Copyright (C) 2013 Rhys Ulerich
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

/** \file
 * \copydoc rand.h
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "rand.h"

// C99 extern declarations for inlined functions from rand.h
extern double tuna_rand_u01(tuna_seed* sd);
extern double tuna_rand_n01(tuna_seed* sd);

tuna_seed tuna_seed_default()
{
    const char * d = getenv("TUNA_SEED");
    unsigned int retval;
    if (!d || 1 != sscanf(d, " %u", &retval)) {
        struct timespec ts;
        clock_gettime(CLOCK_REALTIME, &ts);
        retval = ts.tv_sec + ts.tv_nsec;
    }
    return retval;
}
