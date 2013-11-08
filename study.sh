#!/bin/bash
# Runs a small computational study on algorithmic performance
set -eu
which parallel >/dev/null 2>&1 || (echo 'GNU parallel needed (http://www.gnu.org/software/parallel/)'  && false)
which stats    >/dev/null 2>&1 || (echo 'Program stats needed (https://github.com/rustyrussell/stats)' && false)
parallel --gnu --noswap --nice=10 --eta -j+0                                             \
         --header :                                                                      \
         --results study                                                                 \
         '(for i in $(seq 1 {trials}); do ./examples/basic {invocations}; done) | stats' \
         ::: trials      10 20 50 100 200 500 1000 2000 5000 10000                       \
         ::: invocations 10 20 50 100 200 500 1000 2000 5000 10000
