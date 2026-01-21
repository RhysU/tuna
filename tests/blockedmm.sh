#!/bin/bash
set -eu

echo Default algorithm with 1 thread
export TUNA_ALGO=;             OMP_NUM_THREADS=1 ./examples/blockedmm 128 7
echo ""

echo Degenerate algorithm with 1 thread
export TUNA_ALGO=zero;         OMP_NUM_THREADS=1 ./examples/blockedmm 128 7
echo ""

echo Uniform algorithm with 1 thread
export TUNA_ALGO=uniform;      OMP_NUM_THREADS=1 ./examples/blockedmm 128 7
echo ""

echo Welch1 algorithm with 2 threads
export TUNA_ALGO=welch1;       OMP_NUM_THREADS=2 ./examples/blockedmm 128 7
echo ""

echo Welch1 algorithm, \\nu_ \\to \\infty limit, with 2 threads
export TUNA_ALGO=welch1_nuinf; OMP_NUM_THREADS=2 ./examples/blockedmm 128 7
echo ""

echo Thompson Sampling with 2 threads
export TUNA_ALGO=thompson; OMP_NUM_THREADS=2 ./examples/blockedmm 128 7
echo ""

echo Confidence Bounds with 2 threads
export TUNA_ALGO=ucb1; OMP_NUM_THREADS=2 ./examples/blockedmm 128 7
echo ""
