import matplotlib.pyplot as plt
import pandas as pd
import glob
import math
import os

# Define the tolerance for comparison
tol = 1e-6
return_value = True

# Find all CSV files
csv_files = sorted(glob.glob("FD_check*.csv"))
n_files = len(csv_files)
if n_files == 0:
    raise SystemExit("No files matching 'FD_check*.csv' found.")

# Automatically pick a square-ish grid
ncols = math.ceil(math.sqrt(n_files))
nrows = math.ceil(n_files / ncols)

# Force axes to be a 2D array (even for 1x1), then flatten
fig, axes = plt.subplots(nrows,
                         ncols,
                         figsize=(5 * ncols, 4 * nrows),
                         squeeze=False)
axes = axes.ravel()

for i, csv_file in enumerate(csv_files):
    df = pd.read_csv(csv_file)

    # Main curve
    axes[i].loglog(df["perturbation"].abs(),
                   df["error"].abs(),
                   marker="o",
                   linestyle="-",
                   label="Current")

    # Matching reference: FD_check_foo.csv -> ref_FD_check_foo.csv
    ref_file = csv_file.replace("FD_check", "ref_FD_check", 1)
    ref_file = os.path.join("reference_data", ref_file)
    if os.path.exists(ref_file):
        df_ref = pd.read_csv(ref_file)
        axes[i].loglog(df_ref["perturbation"].abs(),
                       df_ref["error"].abs(),
                       linestyle="--",
                       label="Reference")

        # Check if F and dfdx match the reference
        f = df["F"].values
        dfdx = df["dFdx"].values
        f_ref = df_ref["F"].values
        dfdx_ref = df_ref["dFdx"].values

        rel_error = (f - f_ref) / (f_ref + 1e-30)
        max_rel_error = abs(rel_error).max()
        if max_rel_error > tol:
            print(
                f"Warning: F does not match reference in {csv_file} (max relative error: {max_rel_error})"
            )
            print("Current F:", f)
            print("Reference F:", f_ref)
            return_value = False

        rel_error_dfdx = (dfdx - dfdx_ref) / (dfdx_ref + 1e-30)
        max_rel_error_dfdx = abs(rel_error_dfdx).max()
        if max_rel_error_dfdx > tol:
            print(
                f"Warning: dFdx does not match reference in {csv_file} (max relative error: {max_rel_error_dfdx})"
            )
            print("Current dFdx:", dfdx)
            print("Reference dFdx:", dfdx_ref)
            return_value = False

    # Title from filename
    name = csv_file.split("FD_check_")[-1].replace(".csv", "")
    axes[i].set_title(name)
    axes[i].set_ylabel("Error")
    axes[i].grid(True, which="both", linestyle="--", linewidth=0.5)
    axes[i].legend()

# Hide unused subplots (if any)
for ax in axes[n_files:]:
    ax.axis('off')

# Create plots folder if it does not exist
if not os.path.exists("plots"):
    os.makedirs("plots")

plt.tight_layout()
plt.savefig('plots/FD_comparison.png', dpi=200)
plt.close(fig)

if not return_value:
    raise SystemExit("Discrepancies found compared to reference data.")
