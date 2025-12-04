#!/usr/bin/env python3
import re
import sys
import argparse
import matplotlib.pyplot as plt

def parse_log(path):
    # Parse lines like: "Step =      157 t =   0.1570000E-01"
    # Track whether we're in "Fluid" or "Adjoint fluid" section.
    adjoint_hdr_re = re.compile(r"^\s*-+\s*Adjoint\s+fluid\s*-+\s*$", re.IGNORECASE)
    fluid_hdr_re   = re.compile(r"^\s*-+\s*Fluid\s*-+\s*$", re.IGNORECASE)
    step_time_re   = re.compile(r"Step\s*=\s*(\d+)\s*t\s*=\s*([0-9.+\-Ee]+)")

    mode = "fluid"
    x_idx, tvals, modes, stepnos = [], [], [], []
    count = 0

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if adjoint_hdr_re.search(line):
                mode = "adjoint"
                continue
            if fluid_hdr_re.search(line) and not adjoint_hdr_re.search(line):
                mode = "fluid"
                continue
            m = step_time_re.search(line)
            if m:
                count += 1
                stepno = int(m.group(1))
                t = float(m.group(2).replace("E","e"))
                x_idx.append(count)
                tvals.append(t)
                modes.append(mode)
                stepnos.append(stepno)
    return x_idx, tvals, modes, stepnos

def plot_segments(x, y, modes, out_png, title="Simulation time vs. total timestep index"):
    # Plot as one continuous polyline with color-changing segments per mode.
    # We do this by plotting contiguous runs of the same mode back-to-back.
    plt.figure(figsize=(10,5))
    n = len(x)
    if n == 0:
        raise SystemExit("No step/time entries found.")

    i = 0
    while i < n:
        j = i + 1
        while j < n and modes[j] == modes[i]:
            j += 1
        color = "black" if modes[i] == "fluid" else "red"
        # plot inclusive slice [i, j-1]
        plt.plot(x[i:j], y[i:j], color=color, label=("Forward" if modes[i]=="fluid" else "Adjoint"))
        i = j

    # De-duplicate legend entries
    handles, labels = plt.gca().get_legend_handles_labels()
    uniq = {}
    for h, l in zip(handles, labels):
        uniq[l] = h
    plt.legend(uniq.values(), uniq.keys())

    plt.xlabel("Cumulative timestep")
    plt.ylabel("Simulation time t")
    plt.title(title)
    plt.tight_layout()
    plt.show()
    # plt.savefig(out_png, dpi=150)
    # print(out_png)

def write_csv(x, y, modes, steps, out_csv):
    with open(out_csv, "w", encoding="utf-8") as f:
        f.write("index,step,mode,t\n")
        for i,(xi,si,mi,ti) in enumerate(zip(x,steps,modes,y), start=1):
            f.write(f"{xi},{si},{mi},{ti}\n")
    print(out_csv)

def main():
    ap = argparse.ArgumentParser(description="Plot simulation time vs total steps with segment coloring by mode.")
    ap.add_argument("log", help="Path to log file")
    ap.add_argument("-o", "--out", default="sim_time_vs_total_steps_colored.png", help="Output PNG path")
    ap.add_argument("--csv", default="sim_time_trace.csv", help="Also write a CSV with parsed data")
    args = ap.parse_args()

    x, y, modes, steps = parse_log(args.log)
    write_csv(x, y, modes, steps, args.csv)
    plot_segments(x, y, modes, args.out)

if __name__ == "__main__":
    main()
