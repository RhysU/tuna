#!/bin/bash
set -eu

echo Default algorithm
export TUNA_ALGO=;             ./examples/smallsort
echo ""

echo Welch1 algorithm
export TUNA_ALGO=welch1;       ./examples/smallsort
echo ""

echo Welch1 algorithm, \\nu_ \\to \\infty limit.
export TUNA_ALGO=welch1_nuinf; ./examples/smallsort
echo ""
