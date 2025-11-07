#!/usr/bin/env python3
import re
import argparse
import matplotlib.pyplot as plt

# Regex patterns
rnorm_start_pat = re.compile(
    r"Start step\s+(\d+):\s+rnorm=\s*([+\-]?\d+(?:\.\d*)?(?:[Ee][+\-]?\d+)?)"
)
gmres_init_pat  = re.compile(r"GMRES\(k\)\s+init step")
gmres_inner_pat = re.compile(r"GMRES\(k\)\s+inner step")
gmres_outer_pat = re.compile(r"GMRES\(k\)\s+outer step")
res_pat         = re.compile(r"\|res\|\s*=\s*([+\-]?\d+(?:\.\d*)?(?:[Ee][+\-]?\d+)?)")

def parse_log(path):
    """
    Returns:
      rnorm_steps: list[int]    -- Newton step indices
      rnorm_vals : list[float]  -- rnorm values per Newton step
      gmres_series: list[list[list[float]]]
        Per Newton step -> list of segments; each segment is a list of residuals
        Segment boundaries occur at each GMRES init and after each outer step.
    """
    rnorm_steps, rnorm_vals = [], []
    gmres_series = []         # per-Newton: list of segments; per-segment: list of residuals
    current_segments = None   # list of segments for active Newton
    current_segment  = None   # list of residuals for active segment (resets at init/outer)
    pending_rnorm = []        # queue of (step_idx, rnorm_val) to pair with next GMRES init

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            # Queue rnorm at "Start step N"
            m_rn = rnorm_start_pat.search(line)
            if m_rn:
                step_idx = int(m_rn.group(1))
                rnorm_val = float(m_rn.group(2))
                pending_rnorm.append((step_idx, rnorm_val))
                continue

            # Residual value (used for init/inner/outer)
            m_res = res_pat.search(line)
            if not m_res:
                continue
            res_val = float(m_res.group(1))

            # New Newton's GMRES starts
            if gmres_init_pat.search(line):
                # Flush previous Newton if any
                if current_segments is not None:
                    if current_segment and len(current_segment) > 0:
                        current_segments.append(current_segment)
                    gmres_series.append(current_segments)
                # Start a new Newton
                current_segments = []
                current_segment = []
                # Attach pending rnorm to arrays
                if pending_rnorm:
                    step_idx, rnorm_val = pending_rnorm.pop(0)
                    rnorm_steps.append(step_idx)
                    rnorm_vals.append(rnorm_val)
                # First point of the new segment (local x = 0)
                current_segment.append(res_val)
                continue

            # Inner step: keep adding to current segment
            if gmres_inner_pat.search(line):
                if current_segment is None:
                    current_segment = []
                current_segment.append(res_val)
                continue

            # Outer step: add to current segment, then close segment and start a fresh one
            if gmres_outer_pat.search(line):
                if current_segment is None:
                    current_segment = []
                current_segment.append(res_val)          # include the outer-step residual
                if current_segments is None:
                    current_segments = []
                current_segments.append(current_segment) # close current segment
                current_segment = []                     # new segment starts (local x reset)
                continue

    # Flush trailing data at EOF
    if current_segments is not None:
        if current_segment and len(current_segment) > 0:
            current_segments.append(current_segment)
        gmres_series.append(current_segments)

    return rnorm_steps, rnorm_vals, gmres_series

def main():
    ap = argparse.ArgumentParser(
        description="Two-panel semilogy: rnorm (left), GMRES |res| (right). "
                    "Right panel resets x=0 at GMRES init and after every outer step; color per Newton."
    )
    ap.add_argument("logfile", help="Path to the log file")
    ap.add_argument("-o", "--out", default="lk_twopanels_outer.png",
                    help="Output image filename (default: lk_twopanels_outer.png)")
    ap.add_argument("--title", default="LightKrylov Convergence", help="Figure title")
    args = ap.parse_args()

    rsteps, rvals, gmres_series = parse_log(args.logfile)

    if not rvals:
        raise SystemExit("No rnorm 'Start step' entries found.")
    if not gmres_series:
        raise SystemExit("No GMRES series found (expected 'init step').")

    # Sort rnorm by step index
    pairs = sorted(zip(rsteps, rvals), key=lambda t: t[0])
    xs = [p[0] for p in pairs]
    ys = [p[1] for p in pairs]

    fig, (ax_l, ax_r) = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(args.title)

    # LEFT: rnorm per Newton step
    ax_l.set_yscale("log")
    ax_l.plot(xs, ys, marker="o", linestyle="-", label="rnorm")
    ax_l.set_xlabel("Newton step")
    ax_l.set_ylabel("rnorm")
    ax_l.grid(True, which="both", axis="both", alpha=0.3)
    ax_l.legend(loc="best")

    # RIGHT: GMRES |res|, reset x=0 at init and after each outer step; color per Newton
    ax_r.set_yscale("log")
    color_cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    for i, segments in enumerate(gmres_series, start=1):
        color = color_cycle[(i - 1) % len(color_cycle)]
        first_segment = True
        for seg in segments:
            if not seg:
                continue
            local_x = list(range(len(seg)))  # 0..len-1 for each segment
            label = f"Newton {i}" if first_segment else None
            ax_r.plot(local_x, seg, marker=".", linestyle="-", color=color, alpha=0.95, label=label)
            first_segment = False

    ax_r.set_xlabel("GMRES local iteration (resets at outer steps)")
    ax_r.set_ylabel("|res|")
    ax_r.grid(True, which="both", axis="both", alpha=0.3)
    ax_r.legend(loc="best", title="Series")

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(args.out, dpi=200)
    print(f"Saved plot to {args.out}")

if __name__ == "__main__":
    main()

