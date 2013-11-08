/*
 * Copyright (C) 2013 Rhys Ulerich
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef TUNA_KERNEL_H
#define TUNA_KERNEL_H

/** \file
 * Tracks elapsed time information for compute kernels.
 */

#include <tuna/stats.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Accumulates runtime information about the performance of a compute kernel.
 * Fill storage with zeros, e.g. from POD zero initialization, to construct or
 * reset an instance.
 */
typedef struct tuna_kernel {
    double     outliers[5];  /**< Invariantly-sorted greatest outliers.  */
    tuna_stats stats;        /**< Accumulated statistics sans outliers. */
} tuna_kernel;

/**
 * Record a new elapsed time observation \c t about kernel \c k.
 * If \c t is identically zero, no observation is recorded.
 */
tuna_kernel* tuna_kernel_obs(tuna_kernel * const k, double t);

/**
 * Incorporate all information recorded about kernel \c k into \c s,
 * including any outliers otherwise discarded from consideration.
 */
tuna_stats* tuna_kernel_merge(      tuna_stats  * const s,
                              const tuna_kernel * const k);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* TUNA_KERNEL_H */
