#!/bin/bash

for run in 2382 2421 2463 2648 2711 2827 2834 2868 2898 3000

do

  echo "Running baseline for run ${run}"

  root -l -q "baseline.cc(${run})"

done
