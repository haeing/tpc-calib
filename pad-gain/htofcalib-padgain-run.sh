#!/bin/bash

runs=(2606 2607)

for run in "${runs[@]}"; do
    tmux new-window -t tpc_gain -n "run${run}"

    tmux send-keys -t tpc_gain:"run${run}" \
        "cd ~/work/git/tpc-calib/pad-gain && root -l \"htofcalib-padgain.cc(${run})\"" C-m
done
