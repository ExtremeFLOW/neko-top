#!/usr/bin/env python3
import cmath
import math
import re
import pathlib
import matplotlib.pyplot as plt

# Regex patterns
rnorm_start_pat = re.compile(
    r"Start step\s+(\d+):\s+rnorm=\s*([+\-]?\d+(?:\.\d*)?(?:[Ee][+\-]?\d+)?)"
)
gmres_init_pat  = re.compile(r"GMRES\(k\)\s+init step")
gmres_inner_pat = re.compile(r"GMRES\(k\)\s+inner step")
gmres_outer_pat = re.compile(r"GMRES\(k\)\s+outer step")
res_pat = re.compile(
    r"\|res\|\s*=\s*([+\-]?\d+(?:\.\d*)?(?:[Ee][+\-]?\d+)?)"
)
eig_pat = re.compile(
    r"^\s*(\d+)\s+"
    r"([+\-]?\d+(?:\.\d*)?(?:[Ee][+\-]?\d+)?)\s+"
    r"([+\-]?\d+(?:\.\d*)?(?:[Ee][+\-]?\d+)?)\s+"
    r"([+\-]?\d+(?:\.\d*)?(?:[Ee][+\-]?\d+)?)\s+"
    r"([+\-]?\d+(?:\.\d*)?(?:[Ee][+\-]?\d+)?)\s+([TF])\s*$"
)

DEFAULT_LOG = pathlib.Path("lightkrylov.log")
DEFAULT_EIGS = pathlib.Path("eigs_output.txt")
DEFAULT_CASE = pathlib.Path("newton.case")
DEFAULT_OUT = pathlib.Path("lk_newton_eigs.png")


def resolve_default_path(filename):
    """
    Prefer the current working directory, then fall back to the common
    Newton example locations.
    """
    cwd_path = pathlib.Path.cwd() / filename
    if cwd_path.exists():
        return cwd_path

    repo_root = pathlib.Path(__file__).resolve().parents[2]
    fallback_paths = (
        repo_root / "results" / "newton" / filename,
        repo_root / "logs" / "newton" / filename,
        repo_root / "examples" / "newton" / filename,
    )
    for path in fallback_paths:
        if path.exists():
            return path

    return cwd_path


def parse_end_time(path):
    """
    Extract `case.time.end_time` from the Neko case file.
    """
    pattern = re.compile(
        r'"time"\s*:\s*\{.*?"end_time"\s*:\s*'
        r'([+\-]?\d+(?:\.\d*)?(?:[Ee][+\-]?\d+)?)',
        re.DOTALL,
    )

    text = path.read_text(encoding="utf-8", errors="ignore")
    match = pattern.search(text)
    if not match:
        raise SystemExit(f"Could not find end_time in {path}.")

    return float(match.group(1))

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
    gmres_series = []
    current_segments = None
    current_segment = None
    pending_rnorm = []

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

            # Outer step closes the current segment and resets the local x-axis.
            if gmres_outer_pat.search(line):
                if current_segment is None:
                    current_segment = []
                current_segment.append(res_val)
                if current_segments is None:
                    current_segments = []
                current_segments.append(current_segment)
                current_segment = []
                continue

    # Flush trailing data at EOF
    if current_segments is not None:
        if current_segment and len(current_segment) > 0:
            current_segments.append(current_segment)
        gmres_series.append(current_segments)

    return rnorm_steps, rnorm_vals, gmres_series


def parse_eigs(path):
    """
    Returns the final eigensolve snapshot in `eigs_output.txt`.
    """
    rows = []

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            match = eig_pat.match(line)
            if not match:
                continue
            rows.append(
                {
                    "iter": int(match.group(1)),
                    "re": float(match.group(2)),
                    "im": float(match.group(3)),
                    "modulus": float(match.group(4)),
                    "residual": float(match.group(5)),
                    "converged": match.group(6) == "T",
                }
            )

    if not rows:
        raise SystemExit("No eigenvalue rows found in eigs_output.txt.")

    last_iter = max(row["iter"] for row in rows)
    return [row for row in rows if row["iter"] == last_iter], last_iter


def main():
    log_path = resolve_default_path(DEFAULT_LOG)
    eigs_path = resolve_default_path(DEFAULT_EIGS)
    case_path = resolve_default_path(DEFAULT_CASE)
    out_path = log_path.parent / DEFAULT_OUT

    rsteps, rvals, gmres_series = parse_log(log_path)
    eig_rows, eig_iter = parse_eigs(eigs_path)
    end_time = parse_end_time(case_path)

    if not rvals:
        raise SystemExit("No rnorm 'Start step' entries found.")
    if not gmres_series:
        raise SystemExit("No GMRES series found (expected 'init step').")

    # Sort rnorm by step index
    pairs = sorted(zip(rsteps, rvals), key=lambda t: t[0])
    xs = [p[0] for p in pairs]
    ys = [p[1] for p in pairs]

    fig, (ax_l, ax_r, ax_s) = plt.subplots(1, 3, figsize=(16, 5))
    fig.suptitle("LightKrylov Convergence And Spectrum")

    # LEFT: rnorm per Newton step
    ax_l.set_yscale("log")
    ax_l.plot(xs, ys, marker="o", linestyle="-", label="rnorm")
    ax_l.set_xlabel("Newton step")
    ax_l.set_ylabel("rnorm")
    ax_l.grid(True, which="both", axis="both", alpha=0.3)
    ax_l.legend(loc="best")

    # RIGHT: GMRES |res|, reset x=0 at init and after each outer step.
    ax_r.set_yscale("log")
    color_cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    for i, segments in enumerate(gmres_series, start=1):
        color = color_cycle[(i - 1) % len(color_cycle)]
        first_segment = True
        for seg in segments:
            if not seg:
                continue
            local_x = list(range(len(seg)))
            label = f"Newton {i}" if first_segment else None
            ax_r.plot(
                local_x,
                seg,
                marker=".",
                linestyle="-",
                color=color,
                alpha=0.95,
                label=label,
            )
            first_segment = False

    ax_r.set_xlabel("GMRES local iteration (resets at outer steps)")
    ax_r.set_ylabel("|res|")
    ax_r.grid(True, which="both", axis="both", alpha=0.3)
    ax_r.legend(loc="best", title="Series")

    # SPECTRUM: convert Ritz values to growth rate and frequency.
    spectrum_rows = []
    for row in eig_rows:
        ritz_value = complex(row["re"], row["im"])
        eig_value = cmath.log(ritz_value) / end_time
        spectrum_rows.append(
            {
                "growth_rate": eig_value.real,
                "frequency": eig_value.imag / (2.0 * math.pi),
                "converged": row["converged"],
            }
        )

    conv_rows = [row for row in spectrum_rows if row["converged"]]
    unconverged_rows = [row for row in spectrum_rows if not row["converged"]]

    if conv_rows:
        ax_s.scatter(
            [row["frequency"] for row in conv_rows],
            [row["growth_rate"] for row in conv_rows],
            marker="o",
            label="Converged",
        )
    if unconverged_rows:
        ax_s.scatter(
            [row["frequency"] for row in unconverged_rows],
            [row["growth_rate"] for row in unconverged_rows],
            marker="x",
            label="Unconverged",
        )

    ax_s.axhline(0.0, color="0.7", linewidth=1.0)
    ax_s.axvline(0.0, color="0.7", linewidth=1.0)
    ax_s.set_xlabel("Frequency")
    ax_s.set_ylabel("Growth rate")
    ax_s.set_title(f"Spectrum At Iteration {eig_iter} (T = {end_time:g})")
    ax_s.grid(True, which="both", axis="both", alpha=0.3)
    ax_s.legend(loc="best")

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out_path, dpi=200)
    print(f"Saved plot to {out_path}")

if __name__ == "__main__":
    main()
