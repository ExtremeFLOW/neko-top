#!/bin/bash

# We aim to run on many different numbers of processors. But to avoid edge
# cases, we will only run on prime numbers.
NPROC=(1 2 3 5 7 11 13 17 19 23 29 31 37 41 43 47)

CURRENT_DIR=$(pwd)

TEST_DIR=$(realpath $(dirname $0))
cd $TEST_DIR

genmeshbox 0 1 0 1 0 0.05 20 20 1 .false. .true. .true.

if [ $? -ne 0 ]; then
    echo "genmeshbox failed. Exiting." >&2
    exit 1
fi

# Determine number of physical cores
PROC=$(lscpu | grep "Core(s) per socket" | awk '{print $4}')
echo "Physical processors: $PROC"

for i in "${NPROC[@]}"; do
    [ $i -gt $PROC ] && continue # Skip if we don't have enough processors

    echo "Running with $i processors"
    export NEKO_LOG_FILE=neko_$i.log
    [ -f optimization_data.csv ] && rm -f optimization_data.csv
    mpirun -n $i ./rugby_ball-tests rugby_ball.case

    if [ $? -ne 0 ]; then
        echo "Test failed with $i processors. Exiting." >&2
        exit 1
    fi

    # target=target_$i.csv
    # echo "Comparing output to reference solution $target"
    # if [ -f $target ]; then
    #     diff -q $target optimization_data.csv
    #     if [ $? -ne 0 ]; then
    #         echo "Output $i does not match the reference solution. Exiting." >&2
    #         exit 1
    #     fi
    # fi

    rm -f optimization_data.csv
done

cd $CURRENT_DIR
