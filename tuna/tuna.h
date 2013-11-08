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
 * \mainpage
 *
 * Tuna provides lightweight autotuning over semantically indistinguishable
 * kernels.
 *
 * See the current <a
 * href="https://github.com/RhysU/tuna/blob/master/README.rst">README</a> for a
 * more detailed overview and http://github.com/RhysU/ar for project
 * information.
 */

/** \file
 * TODO
 */

#include <tuna/countof.h>
#include <tuna/kernel.h>
#include <tuna/rand.h>
#include <tuna/stats.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Forward declarations */
struct tuna_state;
typedef struct tuna_state tuna_state;

/**
 * Type signature for all autotuning algorithms.
 *
 * \param[in   ] nkernels How many alternatives are under consideration?
 * \param[inout] kernels  Tracks information about \c nkernels alternatives.
 *                        Must be stored contiguously in memory.
 * \param[inout] seed     Localized pseudo-random number generator state.
 *
 * \return The zero-based index of the kernel that has been selected.
 */
typedef int (*tuna_algorithm)(const int nkernels,
                              tuna_kernel * kernels,
                              tuna_seed * seed);

/** Kernel-independent state required for each autotuning site. */
typedef struct tuna_state {
    tuna_algorithm algo;  /**< The chosen tuning algorithm.   */
    tuna_seed      seed;  /**< Random number generator state. */
} tuna_state;

/**
 * Type signature for all autotuning algorithms.
 *
 * \param[inout] state    Information local to one autotuning site.
 * \param[in   ] nkernels How many alternatives are under consideration?
 * \param[inout] kernels  Tracks information about \c nkernels alternatives.
 *                        Must be stored contiguously in memory.
 *
 * \return The zero-based index of the kernel which should be selected.
 */
static inline int tuna(tuna_state * state,
                       const int nkernels,
                       tuna_kernel * kernels)
{
    return state->algo(nkernels, kernels, &state->seed);
}

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* TUNA_TUNA_H */
