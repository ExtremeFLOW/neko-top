# Static mixer example used for benchmarking

First stage.
- Be able to measure all 3 metrics.
- Then we can define some cases to guide our decisions.

## Weak/strong scaling

We need to make some scaling graphs for both strong and weak scaling. Probably
with finalized memory and IO layout.

## Memory usage

We want to fill out our memory space. So question will be how to balance the IO
with the memory. We want to save enough in memory to retain IO performance, but
as little as possibly to allow a larger case overall.

## IO performance

It is very important to experiment with the number of blocks we work with in the
LUSTRE file system. Striping our files should help performance.