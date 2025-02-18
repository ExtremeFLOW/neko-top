#!/bin/bash

# We aim to run on many different numbers of processors. But to avoid edge
# cases, we will only run on prime numbers.
NPROC=(7 5 3 2 1)

for i in "${NPROC[@]}"; do
    echo "Running with $i processors"
    mpirun -n $i ./neko rugby_ball.case
    if [ -s error.log ]; then
        echo "Error log is not empty. Exiting."
        exit 1
    fi

    target=$RPATH/rugby_verification/target_$i.csv
    echo "Comparing output to reference solution $target"
    if [ -f $target ]; then
        diff -q $target optimization_data.csv
        if [ $? -ne 0 ]; then
            echo "Output $i does not match the reference solution. Exiting." >&2
            exit 1
        fi
    fi

    mv optimization_data.csv target_$i.csv
done
