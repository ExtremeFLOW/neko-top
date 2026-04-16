# Cluster configurations {#clusters}

This section describes the configuration of the Neko-TOP library for different
clusters. Each section will describe the preparation needed to run the library
on the given cluster. The preparation is done in a `prepare.env` file which is 
loaded by the `setup.sh` script. The `prepare.env` file is also loaded by the
`run.sh` script, so it is possible to define environment variables for the
a given example.

## Lumi {#lumi}

LUMI is a supercomputer located in Europe, and it is one of the fastest
supercomputers in the world. The system is based on AMD EPYC processors and
AMD GPUs. The system is managed by the LUMI consortium, which includes
Finland, Sweden, Norway, and Denmark. The system is used for a wide range of
scientific applications, including climate modeling, astrophysics, and
computational fluid dynamics.

Please visit the [LUMI website](https://lumi-supercomputer.eu/) for more
information about the system and its capabilities.

### Limitations

### Environment setup
The following example shows how to set up the environment for LUMI. The
`prepare.env` file should be placed in the root of the Neko-TOP directory. The
`setup.sh` script will automatically load the `prepare.env` file if it exists.
The `prepare.env` file should contain the following lines:

#### CPU Cray
(Tested after LUMI update january 2026)
```bash
# Load modules for the Cray CPU environment
ml CrayEnv cce/20.0.0 cray-mpich/9.0.1 buildtools 2>> /dev/null
export NEKO_CFLAGS="-O3"
export NEKO_FCFLAGS="-O0 -m4"

# # Define the HDF5 support
ml cray-hdf5-parallel/1.14.3.7 2>> /dev/null

# Set CMake build type and compilers
export CMAKE_BUILD_TYPE=Release
export MPIFC=ftn
export MPICC=cc
export MPICXX=CC
export FC=ftn
export CC=cc
```

#### GPU Cray
(Tested after LUMI update january 2026)
```bash
# Load modules for the Cray CPU environment
ml CrayEnv cce/20.0.0 cray-mpich/9.0.1 buildtools 2>> /dev/null
export NEKO_CFLAGS="-O3"
export NEKO_FCFLAGS="-O0 -m4"

# # Define the HDF5 support
ml cray-hdf5-parallel/1.14.3.7 2>> /dev/null

# Load GPU Specific modules
ml craype-accel-amd-gfx90a rocm/6.3.4 2>> /dev/null
export LD_LIBRARY_PATH="$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH"

export HIP_DIR="${ROCM_PATH}"
export HIPCC=hipcc
export NEKO_HIPCC_FLAGS="-O3 --offload-arch=gfx90a"
export NEKO_CONFIG_FLAGS=(--enable-device-mpi)
export MPICH_GPU_SUPPORT_ENABLED=1

# Set CMake build type and compilers
export CMAKE_BUILD_TYPE=Release
export MPIFC=ftn
export MPICC=cc
export MPICXX=CC
export FC=ftn
export CC=cc
```

### Execution of examples
Examples on LUMI can be executed in two different ways, interactively or as a
submitted batch job.

In order to run the examples interactively, you can use the following command:

```bash
salloc --ntasks 2 -t 00:10:00 --partition=small
./run.sh EXAMPLE_NAME
```

Where `EXAMPLE_NAME` is the name of the example you want to run. This command
will allocate 2 tasks (CPI cores) for 10 minutes in the small partition. Then
whenever an example is invoked from the runscript, it will be executed
on the allocated resources.

In order to run the examples as a batch job, you can use the following command:

```bash
./run.sh --submit LUMI EXAMPLE_NAME
```

Where `EXAMPLE_NAME` is the name of the example you want to run. This command
will submit the job to the LUMI batch system based on the jobscript provided in
`scripts/jobscript/lumi/` folder. The jobscript will used to specify options
such as the partition, number of tasks, and time limit. The jobscript will be
submitted to the LUMI batch system using the `sbatch` command.
