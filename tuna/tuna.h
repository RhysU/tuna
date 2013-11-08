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
#include <tuna/ltqnorm.h>
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
 * \param[in   ] nk How many alternatives are under consideration?
 * \param[inout] ks Tracks information about \c nk alternatives.
 *                  Must be stored contiguously in memory.
 * \param[inout] sd Localized pseudo-random number generator state.
 *
 * \return The zero-based index of the kernel that has been selected.
 */
typedef int (*tuna_algorithm)(const int nk,
                              const tuna_kernel * ks,
                              tuna_seed * sd);

/** Kernel-independent state required for each autotuning site. */
typedef struct tuna_state {
    tuna_algorithm al;  /**< The chosen tuning algorithm.   */
    tuna_seed      sd;  /**< Random number generator state. */
} tuna_state;

/**
 * Type signature for all autotuning algorithms.
 *
 * \param[inout] st    Information local to one autotuning site.
 * \param[in   ] nk How many alternatives are under consideration?
 * \param[inout] ks  Tracks information about \c nk alternatives.
 *                        Must be stored contiguously in memory.
 *
 * \return The zero-based index of the kernel which should be selected.
 */
int tuna(tuna_state * st,
         const int nk,
         const tuna_kernel * ks);

/** An autotuning algorithm employing \ref tuna_welch1_nuinf. */
int tuna_algo_welch1_nuinf(const int nk,
                           const tuna_kernel * ks,
                           tuna_seed * sd);

/** An autotuning algorithm employing \ref tuna_welch1. */
int tuna_algo_welch1(const int nk,
                     const tuna_kernel * ks,
                     tuna_seed * sd);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* TUNA_TUNA_H */
