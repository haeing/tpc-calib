#!/bin/bash

runs=(2462 2463 2465 2466)

for run in "${runs[@]}"; do
    tmux new-window -t gain -n "run${run}"

    tmux send-keys -t gain:"run${run}" \
        "cd ~/work/git/tpc-calib/pad-gain && root -l \"physics-padgain.cc(${run})\"" C-m
done
