/**
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

#ifdef __cplusplus
extern "C" {
#endif

/**
 * TODO
 *
 * Must be a POD type.  Likely will contain pointers for runtime dispatch.
 */
typedef struct tuna_algorithm {
} tuna_algorithm;

/**
 * TODO
 *
 * Likely will contain state supporting operations like
 * http://agentzlerich.blogspot.com/2013/03/mergeable-accumulation-of-running-min.html
 *
 * Must be a POD type with zero-initialization providing sane behavior.
 */
typedef struct tuna_kernel {
} tuna_kernel;

/**
 * TODO
 */
typedef struct tuna_sample {
    double start;
    double end;
    int    ndx;
} tuna_sample;

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* TUNA_TUNA_H */
