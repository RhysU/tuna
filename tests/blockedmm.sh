#!/bin/bash
set -eu

echo Default algorithm
export TUNA_ALGO=;             OMP_NUM_THREADS=1 ./examples/blockedmm 64 7
echo ""

echo Degenerate algorithm
export TUNA_ALGO=zero;         OMP_NUM_THREADS=1 ./examples/blockedmm 64 7
echo ""

echo Uniform algorithm
export TUNA_ALGO=uniform;      OMP_NUM_THREADS=1 ./examples/blockedmm 64 7
echo ""

echo Welch1 algorithm
export TUNA_ALGO=welch1;       OMP_NUM_THREADS=2 ./examples/blockedmm 64 7
echo ""

echo Welch1 algorithm, \\nu_ \\to \\infty limit.
export TUNA_ALGO=welch1_nuinf; OMP_NUM_THREADS=2 ./examples/blockedmm 64 7
echo ""
