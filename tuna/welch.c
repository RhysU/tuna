/*
 * Copyright (C) 2013 Rhys Ulerich
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

/** @file
 * @copydoc welch.h
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "welch.h"

// C99 extern declarations for inlined functions from welch.h
extern double tuna_welch1_nuinf(double xA, double sA2, size_t nA,
                                double xB, double sB2, size_t nB);
