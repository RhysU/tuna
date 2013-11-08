/*
 * Copyright (C) 2013 Rhys Ulerich
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef TUNA_RENARD_H
#define TUNA_RENARD_H

/** \file
 * Tabulates the Renard numbers appearing in ISO 3.  See <a
 * href="https://en.wikipedia.org/wiki/Preferred_number"></a> for background.
 */

#ifdef __cplusplus
extern "C" {
#endif

extern const double tuna_R5 [ 5];  /**< Renard's R5  preferred numbers. */
extern const double tuna_R10[10];  /**< Renard's R10 preferred numbers. */
extern const double tuna_R20[20];  /**< Renard's R20 preferred numbers. */
extern const double tuna_R40[40];  /**< Renard's R40 preferred numbers. */
extern const double tuna_R80[80];  /**< Renard's R80 preferred numbers. */

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* TUNA_RENARD_H */
