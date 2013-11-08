/*
 * Copyright (C) 2013 Rhys Ulerich
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

/** \file
 * \copydoc renard.h
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "renard.h"

// Retrieved from https://en.wikipedia.org/wiki/Preferred_number on 8 Nov 2013
const double tuna_R5 [ 5] = {
        1.00, 1.60, 2.50, 4.00, 6.30,
};

// Retrieved from https://en.wikipedia.org/wiki/Preferred_number on 8 Nov 2013
const double tuna_R10[10] = {
        1.00, 1.25, 1.60, 2.00, 2.50,
        3.15, 4.00, 5.00, 6.30, 8.00
};

// Retrieved from https://en.wikipedia.org/wiki/Preferred_number on 8 Nov 2013
const double tuna_R20[20] = {
        1.00, 1.12, 1.25, 1.40, 1.60,
        1.80, 2.00, 2.24, 2.50, 2.80,
        3.15, 3.55, 4.00, 4.50, 5.00,
        5.60, 6.30, 7.10, 8.00, 9.00
};

// Retrieved from https://en.wikipedia.org/wiki/Preferred_number on 8 Nov 2013
const double tuna_R40[40] = {
        1.00, 1.06, 1.12, 1.18, 1.25,
        1.32, 1.40, 1.50, 1.60, 1.70,
        1.80, 1.90, 2.00, 2.12, 2.24,
        2.36, 2.50, 2.65, 2.80, 3.00,
        3.15, 3.35, 3.55, 3.75, 4.00,
        4.25, 4.50, 4.75, 5.00, 5.30,
        5.60, 6.00, 6.30, 6.70, 7.10,
        7.50, 8.00, 8.50, 9.00, 9.50
};

// Retrieved from https://en.wikipedia.org/wiki/Preferred_number on 8 Nov 2013
const double tuna_R80[80] = {
        1.00, 1.03, 1.06, 1.09, 1.12,
        1.15, 1.18, 1.22, 1.25, 1.28,
        1.32, 1.36, 1.40, 1.45, 1.50,
        1.55, 1.60, 1.65, 1.70, 1.75,
        1.80, 1.85, 1.90, 1.95, 2.00,
        2.06, 2.12, 2.18, 2.24, 2.30,
        2.36, 2.43, 2.50, 2.58, 2.65,
        2.72, 2.80, 2.90, 3.00, 3.07,
        3.15, 3.25, 3.35, 3.45, 3.55,
        3.65, 3.75, 3.87, 4.00, 4.12,
        4.25, 4.37, 4.50, 4.62, 4.75,
        4.87, 5.00, 5.15, 5.30, 5.45,
        5.60, 5.80, 6.00, 6.15, 6.30,
        6.50, 6.70, 6.90, 7.10, 7.30,
        7.50, 7.75, 8.00, 8.25, 8.50,
        8.75, 9.00, 9.25, 9.50, 9.75
};
