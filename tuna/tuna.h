/*
 * Copyright (C) 2013 Rhys Ulerich
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef TUNA_TUNA_H
#define TUNA_TUNA_H

/**
 * @mainpage
 *
 * Tuna provides lightweight autotuning over semantically indistinguishable
 * kernels.
 *
 * See the current <a
 * href="https://github.com/RhysU/tuna/blob/master/README.rst">README</a> for a
 * more detailed overview and http://github.com/RhysU/ar for project
 * information.
 */

/** @file
 * Public API
 */

#include <tuna/stats.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Compile-time constants */
enum {

    /** Number of slow outlier samples tracked for each \ref tuna_kernel.  */
    tuna_noutliers = 5

};

/**
 * Accumulates runtime information about the performance of a compute kernel.
 * Fill storage with zeros, e.g. from POD zero initialization, to construct or
 * reset an instance.
 */
typedef struct tuna_kernel {
    double     outliers[tuna_noutliers];
    tuna_stats stats;
} tuna_kernel;

/**
 * Record a new elapsed time observation \c t about kernel \c k.
 * If \c t is identically zero, no observation is recorded.
 */
tuna_kernel* tuna_kernel_obs(tuna_kernel * const k, double t);

/**
 * TODO
 *
 * Must be a POD type.  Likely will contain pointers for runtime dispatch.
 */
typedef struct tuna_algorithm {
} tuna_algorithm;


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* TUNA_TUNA_H */
