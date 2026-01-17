/*
 * Benchmark for tuna_stats operations.
 * Measures performance of tuna_stats_obs and tuna_stats_merge.
 */

#define _POSIX_C_SOURCE 199309L
#include <tuna.h>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/* Prevent compiler from optimizing away results */
static volatile double sink;
static volatile size_t sink_n;

static double
get_time_ns(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec * 1e9 + (double)ts.tv_nsec;
}

/*
 * Benchmark tuna_stats_obs with N observations.
 * Returns nanoseconds per observation.
 */
static double
bench_obs(size_t N)
{
    tuna_stats s = {0};
    double start, end;
    size_t i;

    start = get_time_ns();
    for (i = 0; i < N; ++i) {
        tuna_stats_obs(&s, (double)i * 0.1);
    }
    end = get_time_ns();

    /* Prevent dead code elimination */
    sink = s.m;
    sink_n = s.n;

    return (end - start) / N;
}

/*
 * Benchmark tuna_stats_merge with N merges.
 * Each merge combines a pre-populated src into dst.
 * Returns nanoseconds per merge.
 */
static double
bench_merge(size_t N)
{
    tuna_stats dst = {0};
    tuna_stats src = {0};
    double start, end;
    size_t i;

    /* Prepare src with some observations */
    for (i = 0; i < 1000; ++i) {
        tuna_stats_obs(&src, (double)i * 0.1);
    }
    /* Prepare dst with some observations */
    for (i = 0; i < 1000; ++i) {
        tuna_stats_obs(&dst, (double)i * 0.2);
    }

    start = get_time_ns();
    for (i = 0; i < N; ++i) {
        tuna_stats_merge(&dst, &src);
    }
    end = get_time_ns();

    /* Prevent dead code elimination */
    sink = dst.m;
    sink_n = dst.n;

    return (end - start) / N;
}

/*
 * Benchmark tuna_stats_obs with a near-maximum n value.
 * This exercises the halving path if present.
 * Returns nanoseconds per observation.
 */
static double
bench_obs_near_max(size_t N)
{
    tuna_stats s = {0};
    double start, end;
    size_t i;

    /* Pre-fill to near the halving threshold */
    /* After fix: NMAX = SIZE_MAX/2 - 1, so set n close to it */
    /* For this bench, we simulate by doing fewer obs and checking */
    s.n = (((size_t)-1) / 2) - N - 10;
    s.m = 100.0;
    s.s = 50000.0;

    start = get_time_ns();
    for (i = 0; i < N; ++i) {
        tuna_stats_obs(&s, (double)i * 0.1);
    }
    end = get_time_ns();

    /* Prevent dead code elimination */
    sink = s.m;
    sink_n = s.n;

    return (end - start) / N;
}

/*
 * Benchmark tuna_stats_merge when dst has large n.
 * This exercises the overflow checking path.
 * Returns nanoseconds per merge.
 */
static double
bench_merge_large_n(size_t N)
{
    tuna_stats dst = {0};
    tuna_stats src = {0};
    double start, end;
    size_t i;

    /* Set dst.n to moderate value */
    dst.n = 1000000;
    dst.m = 50.0;
    dst.s = 25000.0;

    /* Set src.n to moderate value */
    src.n = 1000000;
    src.m = 55.0;
    src.s = 30000.0;

    start = get_time_ns();
    for (i = 0; i < N; ++i) {
        tuna_stats_merge(&dst, &src);
    }
    end = get_time_ns();

    /* Prevent dead code elimination */
    sink = dst.m;
    sink_n = dst.n;

    return (end - start) / N;
}

int
main(int argc, char *argv[])
{
    size_t niter_obs = 10000000;
    size_t niter_merge = 1000000;
    size_t niter_near_max = 10000;
    double ns_obs, ns_merge, ns_obs_near_max, ns_merge_large;
    int nruns = 5;
    int run;

    if (argc > 1) {
        nruns = atoi(argv[1]);
    }

    printf("Benchmark: tuna_stats operations\n");
    printf("Iterations: obs=%zu, merge=%zu, near_max=%zu\n\n",
           niter_obs, niter_merge, niter_near_max);

    for (run = 0; run < nruns; ++run) {
        printf("Run %d/%d:\n", run + 1, nruns);

        ns_obs = bench_obs(niter_obs);
        printf("  tuna_stats_obs:        %.3f ns/op\n", ns_obs);

        ns_merge = bench_merge(niter_merge);
        printf("  tuna_stats_merge:      %.3f ns/op\n", ns_merge);

        ns_obs_near_max = bench_obs_near_max(niter_near_max);
        printf("  tuna_stats_obs (max):  %.3f ns/op\n", ns_obs_near_max);

        ns_merge_large = bench_merge_large_n(niter_merge);
        printf("  tuna_stats_merge (lg): %.3f ns/op\n", ns_merge_large);

        printf("\n");
    }

    return 0;
}
