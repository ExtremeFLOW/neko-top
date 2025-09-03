import matplotlib.pyplot as plt
import pandas as pd
import glob
import math

# Find all CSV files
csv_files = sorted(glob.glob("FD_check*.csv"))
n_files = len(csv_files)

# Automatically pick a square-ish grid
ncols = math.ceil(math.sqrt(n_files))
nrows = math.ceil(n_files / ncols)

fig, axes = plt.subplots(nrows, ncols, figsize=(5*ncols, 4*nrows))
axes = axes.flatten()  # makes indexing easier

for i, csv_file in enumerate(csv_files):
    df = pd.read_csv(csv_file)

    axes[i].loglog(df["perturbation"].abs(), df["error"].abs(), marker="o", linestyle="-")

    # Extract title from filename
    name = csv_file.split("FD_check_")[-1].replace(".case.csv", "")
    axes[i].set_title(name)
    axes[i].set_ylabel("Error")
    axes[i].grid(True, which="both", linestyle="--", linewidth=0.5)

# Hide unused subplots
for j in range(i+1, len(axes)):
    axes[j].axis('off')

plt.tight_layout()
plt.savefig('FD_comparison.png')
plt.show()