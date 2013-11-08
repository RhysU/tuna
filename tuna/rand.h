/*
 * Copyright (C) 2013 Rhys Ulerich
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef TUNA_RAND_H
#define TUNA_RAND_H

/** \file
 * Pseudo-random number generation.
 */

#include <stdlib.h>

#include <tuna/ltqnorm.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * What state is required for tuna_u01() and tuna_n01() to be re-entrant safe
 * and orthogonal to all other pseudo-random generators that might be in use?
 */
typedef unsigned int tuna_seed;

/** Generate a uniform draw from <tt>[0, 1]</tt>. */
static inline double tuna_rand_u01(tuna_seed* sd)
{
    return rand_r(sd) / (double) RAND_MAX;
}

/** Generate a draw from <tt>N(0, 1)</tt>. */
static inline double tuna_rand_n01(tuna_seed* sd)
{
    return tuna_ltqnorm(tuna_rand_u01(sd));
}

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* TUNA_RAND_H */
