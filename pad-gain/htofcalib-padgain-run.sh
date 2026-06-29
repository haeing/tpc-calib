#!/bin/bash

runs=(2489 2599 2601 2602 2603 2604 2606 2607)

for run in "${runs[@]}"; do
    tmux new-window -t gain -n "run${run}"

    tmux send-keys -t gain:"run${run}" \
        "cd ~/work/git/tpc-calib/pad-gain && root -l \"htofcalib-padgain.cc(${run})\"" C-m
done
