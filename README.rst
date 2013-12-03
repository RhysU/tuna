Tuna: A Chunk-based Lightweight AutoTuna
========================================

Overview
--------

Ever run into a problem where you can code multiple correct solutions but lack
a robust way to choose the best among them?

Maybe your inputs are small enough that asymptotic complexity isn't the whole
story.  Maybe you lack the time to choose based on careful performance
observations.  Maybe you lack the details to concoct a representative use case
for tuning purposes.  Or maybe you want to do a little future-proofing just in
case the requirements change...

Tuna is for you.  Tuna provides portable, lightweight automatic tuning over
semantically-indistinguishable chunks of code.  It works by misusing the ideas
behind `A/B testing <http://en.wikipedia.org/wiki/A/B_testing>`_ in a
straightforward manner that definitely won't pass peer review.

Example
-------

Hypothetically, say you want to sort many small lists.  Let's autotune over
three candidate ``O(n ln n)`` sorting algorithms with Tuna::

    #include <tuna.h>

    void smallsort(int a[], int n) {
        static tuna_site  si;
        static tuna_chunk ks[3];
        tuna_stack st;
        switch (tuna_pre(&si, &st, ks, tuna_countof(ks))) {
            default: sort_insertion(a, n); break;
            case 1:  sort_qsort    (a, n); break;
            case 2:  sort_heap     (a, n); break;
        }
        tuna_post(&st, ks);
    }

We just wrapped ``sort_insertion()``, ``sort_qsort()`` and ``sort_heap()`` so
that the chunk of logic selected on any invocation of ``smallsort()`` is
dynamically modified based on the observed performance thus far.  Nothing more
is required.

Generally, each call site autotuned by Tuna has a long-lived ``tuna_site`` and
a contiguous array of ``tuna_chunk`` instances, one per chunk under
consideration.  Here, static storage is used to ensure both are
zero-initialized and that they persist across calls.  Additionally, a
``tuna_stack`` shuttles stack-oriented, one-time information from
``tuna_pre()`` to ``tuna_post()``.  Sensible algorithmic defaults are chosen,
but some runtime-selection of behavior can be had.  For details, look in
`tuna.h <tuna/tuna.h>`_ for the ``TUNA_ALGO`` and ``TUNA_SEED`` environment
variables.

This `smallsort example <examples/smallsort.c>`_ is included with Tuna.  Let's
run 1000 sorts on integer lists with 150 elements::

    $ ./examples/smallsort 1000 150
    niter=1000, nelem=150, memory=160 bytes
    insertion              481     2.45536e-05 +/-     3.00301e-06
    qsort(3)                 5     3.90884e-05 +/-     1.18106e-05
    heap                   514     7.19700e-06 +/-     3.18348e-06

The first, second, and third numeric columns are the mean, standard deviation,
and count observed for each chunk, respectively.  Times are given in seconds as
measured by ``CLOCK_PROCESS_CPUTIME_ID``.  On lists of 150 elements,
``sort_insertion()`` is faster and invoked the lion's share of the time we call
``smallsort()``.  The other two chunks are called 5 times each.  Why five?
Tuna omits the three worst outliers from consideration when computing
statistics.  This forgives one-time hiccups like slow startup times.  Two more
calls are required to have a sample standard deviation.

Turning to 165 and then 180 elements per list::

    $ ./examples/smallsort 1000 165
    niter=1000, nelem=165, memory=160 bytes
    insertion              989     1.07414e-05 +/-     7.80852e-06
    qsort(3)                 5     3.97584e-05 +/-     6.18204e-06
    heap                     6     3.14592e-05 +/-     5.43398e-07

    $ ./examples/smallsort 1000 180
    niter=1000, nelem=180, memory=160 bytes
    insertion               41     3.49403e-05 +/-     2.80174e-06
    qsort(3)                 5     4.42012e-05 +/-     4.93605e-06
    heap                   954     1.04687e-05 +/-     6.86004e-06

At 165 elements per list, you can see ``sort_heap()`` has a mean performance
closer to ``sort_insertion()`` now.  Tuna invokes it frequently, sampling it
more frequently than before on account of its standard deviation arguably
making it the fastest given the samples Tuna observed thus far.  At 180
elements, ``sort_heap()`` becomes statistically faster and is therefore more
heavily used.  Tuna discovered change in behavior around 165 samples in the
presence of sampling noise allowing the code to automatically benefit from the
faster chunk regardless of input size.

The algorithmic details are exceedingly simple.  Aside from some basic timing
and running statistics accumulation with outlier tracking, it all boils down to
performing a `one-sided Welch t-test
<http://en.wikipedia.org/wiki/Welch's_t_test>`_ on each call to `tuna_pre()`.
Instead of accepting or rejecting the null hypothesis based on interpreting the
p-value from the t statistic, a uniform random number is drawn on [0,1] to
determine which chunk to select.  This all occurs in ``tuna_algo_welch1()``.
Other algorithms are available by setting ``TUNA_ALGO`` as described by from
reading ``tuna_algo_default()``.  ``tuna_post()`` performs post-invocation
bookkeeping.  Non-time cost measures can be used by calling
``tuna_post_cost()`` instead of ``tuna_post()``.

Possible Use Cases
------------------

Tuna was written to be easy to shoehorn into many similar problem contexts:

1. Should I recompute some mildly expensive value or pay to retrieve it from a
   cache?
2. Should I offload some expensive computation to a coprocessor (Xeon Phi?
   GPU?) or will the offload latency kill me?
3. How many threads should I employ for a compute kernel before resource
   contention causes them to all fall over?
4. Which of several numerics choices will give me the best time-to-solution
   for the particular physics problem I want to solve?
5. Write a decorator for Python to add nice, crisp syntax so you can
   automatically find the fastest of the 57 ways you can write your logic using
   NumPy/SciPy.
6. You tell me.

The necessary ``tuna_site`` and ``tuna_chunk`` data may be stored anywhere.
For example, one could make them member data in a C++ object.  Or protect the
``tuna_pre()`` and ``tuna_post()`` invocations with locking or some thread
local storage.  Or have them live in a large map keyed by some identifier to
permit interrogating what autotuning choices were made by an ``atexit(3)``
hook.

Build and Installation
----------------------

The usual GNU Autotools dance should work::

    ./bootstrap && ./configure --prefix=somewhere && make all check install

Afterwards you can ``include <tuna.h>`` and link with ``-ltuna``.  For those
that hate the GNU Autotools or who simply want to directly incorporate the
functionality, the files `tuna.h <tuna/tuna.h>`_ and `tuna.c <tuna/tuna.c>`_
comprise the entire library and they can be dropped in place nearly anywhere.
