import pandas as pd
import os
import math

# Define the tolerance for comparison
tol = 1e-10
return_value = True

# Define file paths
current_dir = os.path.dirname(os.path.abspath(__file__))
report = open(os.path.join(current_dir, "test_report.txt"), "w")

if not os.path.exists(os.path.join(current_dir, "reference_data")):
    os.makedirs(os.path.join(current_dir, "reference_data"))

for file_name in [
        f for f in os.listdir(current_dir)
        if f.startswith("optimization_data_") and f.endswith(".csv")
]:
    file = os.path.join(current_dir, file_name)

    # Read current and reference data
    df_current = pd.read_csv(file)

    reference_file = os.path.join(current_dir, "reference_data", file_name)

    if not os.path.exists(reference_file):
        df_current.to_csv(reference_file, index=False)
        print(f"Reference data file not found." +
              f"Created new reference file at '{reference_file}'.")
        continue
    df_reference = pd.read_csv(reference_file)

    # Compare data
    for column in df_current.columns:
        if column.split(':')[0].strip() == "backend":
            continue
        if column.split(':')[0].strip() == "subsolver":
            continue

        if column not in df_reference.columns:
            raise SystemExit(f"Column '{column}' not found in reference data.")

    iter = df_current['iter']
    iter_ref = df_reference['iter']
    if not all(iter == iter_ref):
        raise SystemExit(
            "Iteration numbers do not match between current and reference data."
        )

    report.write(f"Comparison Report for {file_name}\n")
    report.write("=" * 80 + "\n\n")

    for column in [
            c for c in df_current.columns
            if c.split(':')[0].strip() != "backend"
            and c.split(':')[0].strip() != "subsolver"
    ]:
        report.write(f"Checking column: {column}\n")
        report.write(f"{'Iter':>6} | {'Current':>15} | {'Reference':>15} | " +
                     f"{'RMSRE':>15} | {'Status':>10}\n")

        for i in range(len(iter)):
            val_current = df_current[column][i]
            val_reference = df_reference[column][i]

            if val_reference == 0.0:
                rmsre = 0.0 if val_current == 0.0 else float("inf")
            elif abs(val_current - val_reference) <= 1e-12:
                rmsre = 0.0
            else:
                rmsre = math.sqrt(
                    ((val_current - val_reference) / val_reference)**2)

            status = "OK" if math.isfinite(rmsre) and rmsre <= tol else "FAIL"
            if status == "FAIL":
                return_value = False

            report.write(
                f"{int(iter[i]):6} | {val_current:15.8e} | " +
                f"{val_reference:15.8e} | {rmsre:15.8e} | {status:>10}\n")
        report.write("\n")

    # Create plots folder if it does not exist
    if not return_value:
        raise SystemExit("Discrepancies found compared to reference data.")

report.close()
