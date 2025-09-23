import matplotlib.pyplot as plt
import pandas as pd
import glob
import math
import os

# Find all CSV files
csv_files = sorted(glob.glob("FD_check*.csv"))
n_files = len(csv_files)
if n_files == 0:
    raise SystemExit("No files matching 'FD_check*.csv' found.")

# Automatically pick a square-ish grid
ncols = math.ceil(math.sqrt(n_files))
nrows = math.ceil(n_files / ncols)

# Force axes to be a 2D array (even for 1x1), then flatten
fig, axes = plt.subplots(nrows, ncols, figsize=(5*ncols, 4*nrows), squeeze=False)
axes = axes.ravel()

for i, csv_file in enumerate(csv_files):
    df = pd.read_csv(csv_file)

    # Main curve
    axes[i].loglog(
        df["perturbation"].abs(),
        df["error"].abs(),
        marker="o",
        linestyle="-",
        label="Current"
    )

    # Matching reference: FD_check_foo.case.csv -> ref_FD_check_foo.case.csv
    ref_file = csv_file.replace("FD_check", "ref_FD_check", 1)
    if os.path.exists(ref_file):
        df_ref = pd.read_csv(ref_file)
        axes[i].loglog(
            df_ref["perturbation"].abs(),
            df_ref["error"].abs(),
            linestyle="--",
            label="Reference"
        )

    # Title from filename
    name = csv_file.split("FD_check_")[-1].replace(".case.csv", "")
    axes[i].set_title(name)
    axes[i].set_ylabel("Error")
    axes[i].grid(True, which="both", linestyle="--", linewidth=0.5)
    axes[i].legend()

# Hide unused subplots (if any)
for ax in axes[n_files:]:
    ax.axis('off')

plt.tight_layout()
plt.savefig('FD_comparison.png', dpi=200)
plt.show()
