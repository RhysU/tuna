#!/bin/bash
set -eu

echo Default algorithm
export TUNA_ALGO=;             ./examples/blockedmm 128 7
echo ""

echo Welch1 algorithm
export TUNA_ALGO=welch1;       ./examples/blockedmm 128 7
echo ""

echo Welch1 algorithm, \\nu_ \\to \\infty limit.
export TUNA_ALGO=welch1_nuinf; ./examples/blockedmm 128 7
echo ""
