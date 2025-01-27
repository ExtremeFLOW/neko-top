#!/bin/bash

# We aim to run on many different numbers of processors. But to avoid edge
# cases, we will only run on prime numbers.
# NPROC=(7 5 3 2 1)
NPROC=(7)

for i in "${NPROC[@]}"; do
    echo "Running with $i processors"
    mpirun -n $i ./neko rugby_ball.case
    if [ -s error.log ]; then
        echo "Error log is not empty. Exiting."
        exit 1
    fi

    if [ ! -f "target_$i.txt" ]; then
        cp optimization_data.txt target_$i.txt
        continue
    fi

    # Check if the run was successful and the output is as expected. Compare
    # the output to the reference solution.
    diff -q target_$i.txt optimization_data.txt
    if [ $? -ne 0 ]; then
        echo "Output does not match the reference solution. Exiting." >&2
        exit 1
    fi

done
