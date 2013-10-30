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
 * TODO
 */

#include <tuna/kernel.h>
#include <tuna/ltqnorm.h>
#include <tuna/stats.h>

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


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* TUNA_TUNA_H */
