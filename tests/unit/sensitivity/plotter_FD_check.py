import pandas as pd
import matplotlib.pyplot as plt
import glob
import os

# Find all CSV files starting with "FD_check"
csv_files = sorted(glob.glob("FD_check*.case.csv"))

if not csv_files:
    raise FileNotFoundError("No files starting with 'FD_check' found in the current directory.")

# Prepare figure with one subplot per file
fig, axes = plt.subplots(len(csv_files), 1, figsize=(6, 4 * len(csv_files)), constrained_layout=True)

# Ensure axes is iterable
if len(csv_files) == 1:
    axes = [axes]

for ax, file in zip(axes, csv_files):
    # Extract problem name from file
    base = os.path.basename(file)
    problem_name = base[len("FD_check_") : -len(".case.csv")]
    
    # Load CSV, skip the first entry
    df = pd.read_csv(file).iloc[1:]
    
    # Plot absolute error in log-log scale
    ax.loglog(df["perturbation"], df["error"].abs(), marker="o", linestyle="-")
    
    # Labeling
    ax.set_xlabel("Perturbation")
    ax.set_ylabel("|Error|")
    ax.set_title(problem_name)
    ax.grid(True, which="both", linestyle="--", linewidth=0.5)

plt.show()

