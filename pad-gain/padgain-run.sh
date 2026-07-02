#!/bin/bash

runs=(2450 2451 2452 2453 2454 2456 2457)

for run in "${runs[@]}"; do
    tmux new-window -t gain -n "run${run}"

    tmux send-keys -t gain:"run${run}" \
        "cd ~/work/git/tpc-calib/pad-gain && root -l \"physics-padgain.cc(${run})\"" C-m
done
