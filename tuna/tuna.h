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

#include <time.h>

#include <tuna/algo.h>
#include <tuna/countof.h>
#include <tuna/kernel.h>
#include <tuna/ltqnorm.h>
#include <tuna/rand.h>
#include <tuna/stats.h>

#ifdef __cplusplus
extern "C" {
#endif

/** Kernel-independent state required for each autotuning site. */
typedef struct tuna_state {
    tuna_algo       al;  /**< The chosen tuning algorithm.                */
    tuna_seed       sd;  /**< Random number generator state.              */
    int             ik;  /**< Index of the most recently selected kernel. */
    struct timespec ts;  /**< Records clock_gettime(2) in tuna_pre().     */
} tuna_state;

/**
 * Invoke the currently selected autotuning algorithm.
 *
 * \param[inout] st Information local to one autotuning site.
 * \param[in   ] ks Tracks information about \c nk alternatives.
 *                  Must be stored contiguously in memory.
 * \param[in   ] nk How many alternatives are under consideration?
 *
 * \return The zero-based index of the kernel which should be selected.
 */
int tuna_pre(tuna_state* st,
             const tuna_kernel* ks,
             const int nk);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* TUNA_TUNA_H */
