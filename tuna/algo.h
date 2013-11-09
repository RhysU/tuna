/*
 * Copyright (C) 2013 Rhys Ulerich
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef TUNA_ALGO_H
#define TUNA_ALGO_H

/** \file
 * Available autotuning algorithms.
 */

#include <tuna/kernel.h>
#include <tuna/rand.h>

#ifdef __cplusplus
extern "C" {
#endif

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
typedef int (*tuna_algo)(const int nk,
                         const tuna_kernel* ks,
                         tuna_seed* sd);

/** An autotuning algorithm employing \ref tuna_welch1_nuinf. */
int tuna_algo_welch1_nuinf(const int nk,
                           const tuna_kernel* ks,
                           tuna_seed* sd);

/** An autotuning algorithm employing \ref tuna_welch1. */
int tuna_algo_welch1(const int nk,
                     const tuna_kernel* ks,
                     tuna_seed* sd);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* TUNA_ALGO_H */
