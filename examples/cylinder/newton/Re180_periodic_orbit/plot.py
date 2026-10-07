
#!/usr/bin/env python3
import re
import pathlib
import matplotlib.pyplot as plt

# Regex patterns
rnorm_start_pat = re.compile(
    r"Start step\s+(\d+):\s+rnorm=\s*([+\-]?\d+(?:\.\d*)?(?:[Ee][+\-]?\d+)?)"
)
gmres_init_pat = re.compile(r"GMRES\(k\)\s+init step")
gmres_inner_pat = re.compile(r"GMRES\(k\)\s+inner step")
gmres_outer_pat = re.compile(r"GMRES\(k\)\s+outer step")
res_pat = re.compile(
    r"\|res\|\s*=\s*([+\-]?\d+(?:\.\d*)?(?:[Ee][+\-]?\d+)?)"
)

DEFAULT_LOG = pathlib.Path("lightkrylov.log")
DEFAULT_OUT = pathlib.Path("lk_newton_convergence.png")


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


def parse_log(path):
    """
    Returns:
      rnorm_steps: list[int]
      rnorm_vals: list[float]
      gmres_series: list[list[list[float]]]

    GMRES segments reset at each init and after each outer step.
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

            # Residual value
            m_res = res_pat.search(line)
            if not m_res:
                continue
            res_val = float(m_res.group(1))

            # New Newton GMRES starts
            if gmres_init_pat.search(line):
                if current_segments is not None:
                    if current_segment:
                        current_segments.append(current_segment)
                    gmres_series.append(current_segments)

                current_segments = []
                current_segment = []

                if pending_rnorm:
                    step_idx, rnorm_val = pending_rnorm.pop(0)
                    rnorm_steps.append(step_idx)
                    rnorm_vals.append(rnorm_val)

                current_segment.append(res_val)
                continue

            # Inner step
            if gmres_inner_pat.search(line):
                if current_segment is None:
                    current_segment = []
                current_segment.append(res_val)
                continue

            # Outer step closes the current segment
            if gmres_outer_pat.search(line):
                if current_segment is None:
                    current_segment = []
                current_segment.append(res_val)

                if current_segments is None:
                    current_segments = []

                current_segments.append(current_segment)
                current_segment = []
                continue

    # Flush trailing data
    if current_segments is not None:
        if current_segment:
            current_segments.append(current_segment)
        gmres_series.append(current_segments)

    return rnorm_steps, rnorm_vals, gmres_series


def main():
    log_path = resolve_default_path(DEFAULT_LOG)
    out_path = log_path.parent / DEFAULT_OUT

    rsteps, rvals, gmres_series = parse_log(log_path)

    if not rvals:
        raise SystemExit("No rnorm 'Start step' entries found.")
    if not gmres_series:
        raise SystemExit("No GMRES series found (expected 'init step').")

    # Sort rnorm by Newton step index
    pairs = sorted(zip(rsteps, rvals), key=lambda t: t[0])
    xs = [p[0] for p in pairs]
    ys = [p[1] for p in pairs]

    # Two-panel figure
    fig, (ax_l, ax_r) = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle("LightKrylov Newton and GMRES Convergence")

    # LEFT: Newton convergence
    ax_l.set_yscale("log")
    ax_l.plot(
        xs, ys,
        marker="o",
        linestyle="-",
        label="rnorm"
    )
    ax_l.set_xlabel("Newton step")
    ax_l.set_ylabel("rnorm")
    ax_l.grid(True, which="both", axis="both", alpha=0.3)
    ax_l.legend(loc="best")

    # RIGHT: GMRES convergence
    # Local iteration resets at init and after each outer step.
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

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out_path, dpi=200)
    print(f"Saved plot to {out_path}")


if __name__ == "__main__":
    main()

