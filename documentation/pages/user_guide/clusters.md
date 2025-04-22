# Cluster configurations {#clusters}


## Lumi-C GNU Compiler{#lumi}
prepare.env

```bash
module purge 2>/dev/null
ml LUMI/24.03 PrgEnv-gnu cray-libsci buildtools

export MPICC=$(which cc)
export MPICXX=$(which CC)
export MPIFC=$(which ftn)
export FC=$(which ftn)
export CXX=$(which CC)
export CC=$(which cc)
```