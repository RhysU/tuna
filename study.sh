#!/bin/bash
# Runs a small computational study on algorithmic performance
set -e
which parallel >/dev/null 2>&1 || (echo 'GNU parallel needed (http://www.gnu.org/software/parallel/)'  && false)
which stats    >/dev/null 2>&1 || (echo 'Program stats needed (https://github.com/rustyrussell/stats)' && false)
ntrial=$1; echo -n "Running ${ntrial:=1e3} trials"
niter=$2;  echo    " each for {${niter:=1e4 5e3 2e3 1e3 5e2 2e2 1e2 5e1 2e1 1e1}} iterations"
parallel --gnu --noswap --nice=10 --eta -j+0                                               \
         '(for i in $(seq 1 {1}); do ./examples/basic {2}; done) | stats > result.{1}.{2}' \
                  ::: $ntrial ::: $niter
