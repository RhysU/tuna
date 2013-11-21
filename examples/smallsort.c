/**
 * Copyright (C) 2013 Rhys Ulerich
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

/** \file
 * A Tuna example for dynamically investigating sorting tradeoffs for small
 * lists.  This provides an autotuned way to account for of compiler- and
 * optimization-dependent behavior.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <tuna.h>

static const char *names[] = { "insertion", "qsort(3)", "heap" };
void sort_insertion(int a[], int array_size);
void sort_qsort    (int a[], int array_size);
void sort_heap     (int a[], int array_size);

static tuna_site  si;     // Normally si and ks would be inside smallsort()
static tuna_chunk ks[3];  // but they are global to permit querying them
void smallsort(int a[], int n) {
    tuna_stack st;
    switch (tuna_pre(&si, &st, ks, tuna_countof(ks))) {
        default: sort_insertion(a, n); break;
        case 1:  sort_qsort    (a, n); break;
        case 2:  sort_heap     (a, n); break;
    }
    tuna_post(&st, ks);
}

int main(int argc, char *argv[])
{
    struct timespec tp;
    clock_gettime(CLOCK_REALTIME, &tp);
    srand((unsigned int) (tp.tv_sec + tp.tv_nsec));

    // Parse any incoming command line arguments
    const int niter = argc > 1 ? atof(argv[1]) : 1000; // Iteration count?
    const int nelem = argc > 2 ? atof(argv[2]) :   10; // Elements to sort?

    int data[nelem];
    for (int i = 0; i < niter; ++i) {

        for (int j = 0; j < nelem; ++j) { // Random input
            data[j] = rand();
        }
        smallsort(data, nelem);           // Autotuned
    }

    // Display settings and static memory overhead required for autotuning
    printf("niter=%d, nelem=%d, memory=%zd bytes\n",
           niter, nelem, sizeof(si) + sizeof(ks));

    // Display observations from each alternative
    for (int i = 0; i < tuna_countof(ks); ++i) {
        tuna_stats o = {};
        tuna_chunk_merge(&o, ks + i);
        printf("%12s\tm=%g,\ts=%g,\tc=%zd\n", names[i],
            tuna_stats_avg(&o), tuna_stats_std(&o), tuna_stats_cnt(&o));
    }

    return EXIT_SUCCESS;
}

// From http://www.codebeach.com/2008/09/sorting-algorithms-in-c.html
void sort_insertion(int a[], int array_size)
{
     int i, j, index;
     for (i = 1; i < array_size; ++i)
     {
          index = a[i];
          for (j = i; j > 0 && a[j-1] > index; j--)
               a[j] = a[j-1];

          a[j] = index;
     }
}

// From http://www.codebeach.com/2008/09/sorting-algorithms-in-c.html
void down_heap(int a[], int root, int bottom)
{
     int maxchild, temp, child;
     while (root*2 < bottom)
     {
          child = root * 2 + 1;
          if (child == bottom)
          {
               maxchild = child;
          }
          else
          {
               if (a[child] > a[child + 1])
                    maxchild = child;
               else
                    maxchild = child + 1;
          }

          if (a[root] < a[maxchild])
          {
               temp = a[root];
               a[root] = a[maxchild];
               a[maxchild] = temp;
          }
          else return;

          root = maxchild;
     }
}

// From http://www.codebeach.com/2008/09/sorting-algorithms-in-c.html
void sort_heap(int a[], int array_size)
{
     int i;
     for (i = (array_size/2 -1); i >= 0; --i)
     {
          down_heap(a, i, array_size-1);
     }

     for (i = array_size-1; i >= 0; --i)
     {
          int temp;
          temp = a[i];
          a[i] = a[0];
          a[0] = temp;
          down_heap(a, 0, i-1);
     }
}

int qsort_compar(const void *a, const void *b)
{
    return *(const int*)a - *(const int*)b;
}

void sort_qsort(int a[], int array_size)
{
    return qsort(a, array_size, sizeof(int), &qsort_compar);
}
