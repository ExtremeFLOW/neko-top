#%%
# Import the necessary libraries and set up the environment

import os
import sys
import glob
import numpy as np
import re
import matplotlib.pyplot as plt

# Define some temporary variables which should be replaced by input parameters
scaling_example = "scaling_rugby"
title = "Rugby ball"

# Define globally useful variables
MAIN_DIR = os.path.dirname(
    os.path.abspath(os.path.dirname(os.path.abspath(__file__)) + "/.."))
example_directory = MAIN_DIR + "/examples"
result_directory = MAIN_DIR + "/results/lumi"

#%%
# Get the list of directories in the example directory
example_dir = os.path.join(example_directory, scaling_example)
result_dir = os.path.join(result_directory, scaling_example)
if not os.path.exists(example_dir):
    print(f"Example directory {example_dir} does not exist.")
    sys.exit(1)
if not os.path.exists(result_dir):
    print(f"Result directory {result_dir} does not exist.")
    sys.exit(1)

# List the case files in the example directory
cases = [
    d for d in glob.glob(os.path.join(example_dir, "*"))
    if os.path.isfile(d) and d.endswith(".case")
]

# Extract the number of processes from the basename of the case files
nprocs = [int(os.path.basename(d).split(".")[0]) for d in cases]

# Sort the cases based on the number of processors
nprocs, cases = zip(*sorted(zip(nprocs, cases)))
nprocs = np.array(nprocs, dtype=int)

#%%
# Initialize the arrays to store the time spent in each experiment

N = len(cases)
TotalTime = -1.0 * np.ones(N)
FluidTime = -1.0 * np.ones(N)
AdjointTime = -1.0 * np.ones(N)
MMATime = -1.0 * np.ones(N)

ex = "(\\d+\\.\\d+E[+-]\\d+)"
# Define the regex patterns to extract the time spent in each experiment
re_total = re.compile(f"Execution time:\\s*(\\d+):(\\d+):(\\d+)")
re_fluid = re.compile(f"Fluid\\s*({ ex })\\s*({ ex })\\s*({ ex })")
re_adjoint = re.compile(f"Adjoint\\s*({ ex })\\s*({ ex })\\s*({ ex })")
re_mma = re.compile(f"MMA\\s*({ ex })\\s*({ ex })\\s*({ ex })")

# Go through each experiment and extract the time spent
for i in range(N):
    # Get the case file
    case_file = cases[i]
    experiment_name = os.path.basename(case_file).split(".")[0]

    if not os.path.exists(os.path.join(result_dir, experiment_name)):
        continue

    # Use regex to extract the time spent from the output.log file
    output_log = os.path.join(result_dir, experiment_name, "output.log")
    neko_log = os.path.join(result_dir, experiment_name,
                            experiment_name + ".log")

    if not os.path.exists(output_log):
        print(f"Log file {output_log} does not exist.")
        sys.exit(1)

    # Use regex to find the total_pattern "Execution time: hh:mm:ss" in the
    # output.log file
    with open(output_log, "r") as f:
        content = f.read()
        match = re_total.search(content)
        if match:
            # Extract hours, minutes, and seconds from the regex groups
            h, m, s = map(int, match.groups())
            # Convert the time to seconds
            TotalTime[i] = h * 3600.0 + m * 60.0 + s

    # Use regex to parse the time spent in the fluid solver from the neko log
    with open(neko_log, "r") as f:
        content = f.read()

        # Extract the time spent
        fluid_match = re_fluid.search(content)
        adjoint_match = re_adjoint.search(content)
        mma_match = re_mma.search(content)

        # Save the time spent in each experiment
        if fluid_match: FluidTime[i] = float(fluid_match.group(1))
        if adjoint_match: AdjointTime[i] = float(adjoint_match.group(1))
        if mma_match: MMATime[i] = float(mma_match.group(1))

#%%
# Construct the plots

# Remove the cases where the time is -1.0
Valid = np.where(TotalTime > 0.0)[0]
N = len(Valid)
nprocs = nprocs[Valid]
TotalTime = TotalTime[Valid]
FluidTime = FluidTime[Valid]
AdjointTime = AdjointTime[Valid]
MMATime = MMATime[Valid]

# Report the time spent in each experiment
print("Time spent in each experiment:")
for i in range(N):
    print(f"  {nprocs[i]} procs: {TotalTime[i]:.2f} s")
    if FluidTime[i] > 0.0:
        print(f"    Fluid time: {FluidTime[i]:.2f} s")
    if AdjointTime[i] > 0.0:
        print(f"    Adjoint time: {AdjointTime[i]:.2f} s")
    if MMATime[i] > 0.0:
        print(f"    MMA time: {MMATime[i]:.2f} s")
    print(f"  Total time: {TotalTime[i]:.2f} s")

# Plot the total time spent in each experiment
plt.figure(figsize=(5, 5))
# Plot the total time spent in each experiment
plt.plot(nprocs,
         TotalTime[0] / (nprocs / nprocs[0]),
         "--",
         label="Ideal",
         linewidth=0.75)
plt.plot(nprocs, TotalTime, "o-", label="Total to solution")
if np.any(FluidTime > 0.0):
    plt.plot(nprocs, FluidTime, "o-", label="Forward time")
if np.any(AdjointTime > 0.0):
    plt.plot(nprocs, AdjointTime, "o-", label="Sensitivity time")
if np.any(MMATime > 0.0):
    plt.plot(nprocs, MMATime, "o-", label="MMA time")

# Set a title and labels
plt.title(f"Strong scaling of {title} Optimization")

plt.xlabel("Number of processes")
plt.xscale("log", base=2)
plt.xticks(nprocs, nprocs)

plt.ylabel("Time (s)")

# Set the y limits to be rounded to the nearest 10^n
ymin = np.min(TotalTime[0] / (nprocs / nprocs[0]))
ymax = np.max(TotalTime)
ymin = 10**np.floor(np.log10(ymin))
ymax = 10**np.ceil(np.log10(ymax))
plt.ylim(ymin, ymax)

plt.yscale("log", base=10)
# plt.yticks(nprocs, nprocs)

plt.legend()
plt.grid()

# Save the plot
plt.savefig(os.path.join(result_dir, "strong_scaling.png"))
