#!/usr/bin/env bash
set -eo pipefail

# Run the cylinder examples in dependency order without changing the checked-in
# case files. Every stage is executed through scripts/functions.sh from inside
# examples/cylinder/gauntlet.

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
CYLINDER="$ROOT/examples/cylinder"
GAUNTLET="$CYLINDER/gauntlet"
SHARED="$GAUNTLET/shared"
WORK="$GAUNTLET/work"
RESULTS="$GAUNTLET/results"
RUN_3D=${RUN_3D:-true}
NPROCS=${NPROCS:-1}
SOLVER_TOL=${GAUNTLET_SOLVER_TOL:-1.0e-7}
ANALYSIS_TOL=${GAUNTLET_ANALYSIS_TOL:-1.0e-8}
NEWTON_TOL=${GAUNTLET_NEWTON_TOL:-1.0e-6}
SPAN_LENGTH=${GAUNTLET_SPAN_LENGTH:-4}
SPAN_LAYERS=${GAUNTLET_SPAN_LAYERS:-4}

usage() {
  cat <<EOF
Usage: $0 [--clean] [--no-3d] [--span-length LENGTH] [--span-layers N]

Cases with a non-empty directory under gauntlet/results are skipped. Use
--clean to remove the gauntlet directory and rerun every selected case.
Changing the 3D spanwise parameters requires --clean if previous 3D
results exist, to prevent mixing results from incompatible meshes.

3D geometry options:
  --span-length LENGTH    spanwise domain length (default: 4)
  --span-layers N         number of spanwise elements (default: 4)

Environment overrides:
  NPROCS                  MPI ranks per case (default: 1)
  GAUNTLET_SOLVER_TOL    fluid absolute tolerance (default: 1e-7)
  GAUNTLET_ANALYSIS_TOL  eigensolver/SVD tolerance (default: 1e-8)
  GAUNTLET_NEWTON_TOL    Newton/periodic-orbit tolerance (default: 1e-6)
  RUN_3D=false            skip the 3D Floquet cases
  GAUNTLET_SPAN_LENGTH    default for --span-length (default: 4)
  GAUNTLET_SPAN_LAYERS    default for --span-layers (default: 4)
EOF
}

CLEAN=false
while [[ $# -gt 0 ]]; do
  case "$1" in
    --clean) CLEAN=true; shift ;;
    --no-3d) RUN_3D=false; shift ;;
    --span-length|--length)
      [[ $# -ge 2 ]] || { echo "Missing value for $1" >&2; exit 2; }
      SPAN_LENGTH=$2; shift 2 ;;
    --span-layers|--layers)
      [[ $# -ge 2 ]] || { echo "Missing value for $1" >&2; exit 2; }
      SPAN_LAYERS=$2; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

# Validate and canonicalize the values (e.g. 4 and 4.0 mean the same thing).
validated_span=$(python3 - "$SPAN_LENGTH" "$SPAN_LAYERS" <<'PY_VALIDATE'
import math
import sys

try:
    length = float(sys.argv[1])
    layers = int(sys.argv[2])
    if not math.isfinite(length) or length <= 0 or layers <= 0:
        raise ValueError("spanwise length and layer count must be positive")
except ValueError as exc:
    print(f"Invalid spanwise parameters: {exc}", file=sys.stderr)
    sys.exit(2)

print(f"{length:.15g} {layers}")
PY_VALIDATE
) || exit 2
read -r SPAN_LENGTH SPAN_LAYERS <<< "$validated_span"

# Initialize the project environment exactly as run.sh does. This must happen
# before mesh preparation, since the mesh converters need Neko's dependencies.
initialize_environment() {
  export CURRENT_DIR=$(pwd)
  export MAIN_DIR="$ROOT"
  export EXTERNAL_DIR="$MAIN_DIR/external"

  if [[ -f "$MAIN_DIR/prepare.env" ]]; then
    source "$MAIN_DIR/prepare.env"
  fi

  export EPATH="$MAIN_DIR/examples"
  export RPATH="$MAIN_DIR/results"
  export LPATH="$MAIN_DIR/logs"
  export SPATH="$MAIN_DIR/scripts/"
  export DPATH="$MAIN_DIR/data"
  export DLPATH="$MAIN_DIR/data_local"

  if [[ -z "${NEKO_DIR:-}" ]]; then
    export NEKO_DIR="$MAIN_DIR/external/neko"
  fi
  export NEKO_DIR=$(realpath "$NEKO_DIR")

  source "$MAIN_DIR/scripts/dependencies.sh"
  find_json_fortran
  find_neko
}

initialize_environment

if [[ "$CLEAN" == true ]]; then
  rm -rf "$GAUNTLET"
fi

mkdir -p "$SHARED" "$WORK"

# Results from another spanwise discretization are not interchangeable.
# Older versions of this script hard-coded length=4 and layers=4, so their
# results can be identified even without a configuration stamp.
SPAN_CONFIG_FILE="$SHARED/3d_spanwise_config.txt"
if [[ "$RUN_3D" == true ]]; then
  requested_config="$SPAN_LENGTH $SPAN_LAYERS"
  previous_config=
  if [[ -f "$SPAN_CONFIG_FILE" ]]; then
    previous_config=$(cat "$SPAN_CONFIG_FILE")
  elif [[ -d "$RESULTS/Re180/LNS_floquete_3D" || \
          -d "$RESULTS/Re180/adjoint_floquete_3D" || \
          -d "$RESULTS/Re200/LNS_floquete_3D" || \
          -d "$RESULTS/Re200/adjoint_floquete_3D" ]]; then
    previous_config="4 4"
  fi
  if [[ -n "$previous_config" && "$requested_config" != "$previous_config" ]]; then
    echo "3D results/mesh were prepared with span parameters: $previous_config" >&2
    echo "Requested parameters: $requested_config" >&2
    echo "Use --clean to rerun without mixing incompatible 3D results." >&2
    exit 2
  fi
fi

prepare_2d_mesh() {
  if [[ ! -f "$SHARED/1cyl.nmsh" ]]; then
    cp "$CYLINDER/Re40/baseflow/newton.geo" "$SHARED/newton.geo"
    cat > "$SHARED/prepare.sh" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
gmsh -0 newton.geo >/dev/null
gmsh2nek <<'MESH_EOF' >/dev/null
2
newton
0
MESH_EOF
rea2nbin newton.re2 1cyl.nmsh
EOF
    chmod u+x "$SHARED/prepare.sh"
    cp "$ROOT/scripts/functions.sh" "$SHARED/functions.sh"
    (
      cd "$SHARED"
      export MAIN_DIR="$ROOT"
      export NEKO_DIR="${NEKO_DIR:-$ROOT/external/neko}"
      source functions.sh
      prepare
    )
  fi
}

patch_case() {
  local input=$1
  local output=$2
  local case_label=$3

  jq --arg solver "$SOLVER_TOL" \
     --arg analysis "$ANALYSIS_TOL" \
     --arg newton "$NEWTON_TOL" \
     --arg case_label "$case_label" '
    def number($x): ($x | tonumber);
    .case.mesh_file = "newton.nmsh"
    | walk(if type == "object" then
      if has("absolute_tolerance") then
        .absolute_tolerance = number($solver)
      elif has("relative_tolerance") then
        .relative_tolerance = number($newton)
      elif has("tolerance") then
        .tolerance = number($analysis)
      else . end
    else . end)
    | if .case.lightkrylov.newton? then
        .case.lightkrylov.newton.relative_tolerance = number($newton)
      else . end
    | if .case.lightkrylov.periodic_orbit? then
        .case.lightkrylov.periodic_orbit.absolute_tolerance = number($newton)
        | .case.lightkrylov.periodic_orbit.outer_iteration_output_prefix =
            ("outer_periodic_orbit_" + $case_label)
        | .case.lightkrylov.periodic_orbit.movie_output_prefix =
            ("periodic_orbit_" + $case_label + "_movie")
      else . end
  ' "$input" > "$output"
}

case_has_results() {
  local results_dir="$RESULTS/$1"

  [[ "$CLEAN" != true ]] && [[ -d "$results_dir" ]] || return 1

  case "$1" in
    Re180/periodic_orbit)
      [[ -f "$results_dir/baseflow/baseflow0.f00000" \
        && -f "$results_dir/baseflow/baseflow0.nek5000" ]]
      ;;
    Re180/LNS_floquete)
      [[ -f "$results_dir/floquet_spectrum.csv" ]]
      ;;
    Re180/adjoint_floquete)
      [[ -f "$results_dir/adjoint_floquet_spectrum.csv" ]]
      ;;
    Re180/LNS_floquete_3D)
      [[ -f "$results_dir/floquet_spectrum.csv" ]]
      ;;
    Re180/adjoint_floquete_3D)
      [[ -f "$results_dir/adjoint_floquet_3D_spectrum.csv" ]]
      ;;
    Re200/periodic_orbit)
      [[ -f "$results_dir/baseflow/baseflow0.f00000" \
        && -f "$results_dir/baseflow/baseflow0.nek5000" ]]
      ;;
    Re200/LNS_floquete_3D)
      [[ -f "$results_dir/floquet_spectrum.csv" ]]
      ;;
    Re200/adjoint_floquete_3D)
      [[ -f "$results_dir/adjoint_floquet_3D_spectrum.csv" ]]
      ;;
    *)
      find "$results_dir" -mindepth 1 -print -quit | grep -q .
      ;;
  esac
}

stage_case() {
  local name=$1
  local source_dir=$2
  local source_case=$3
  local binary=$4
  local mesh=${5:-2d}
  local dir="$WORK/$name"

  case_has_results "$name" && return 0

  rm -rf "$dir"
  mkdir -p "$dir"
  patch_case "$source_dir/$source_case" "$dir/$source_case" "${name%%/*}"
  cp "$binary" "$dir/neko"
  cp "$ROOT/scripts/functions.sh" "$dir/functions.sh"
  chmod u+x "$dir/neko"

  if [[ "$mesh" == 2d ]]; then
    cat > "$dir/prepare.sh" <<EOF
#!/usr/bin/env bash
set -euo pipefail
cp "$SHARED/1cyl.nmsh" newton.nmsh
EOF
  else
    cat > "$dir/prepare.sh" <<EOF
#!/usr/bin/env bash
set -euo pipefail
cp "$SHARED/1cyl_3d.nmsh" newton.nmsh
EOF
  fi
  chmod u+x "$dir/prepare.sh"
}

run_staged_case() {
  local name=$1
  local source_case=$2
  local dir="$WORK/$name"

  if case_has_results "$name"; then
    echo "[gauntlet] skipping $name (results already exist)"
    return 0
  fi

  echo "[gauntlet] $name"
  (
    cd "$dir"
    export MAIN_DIR="$ROOT"
    export RPATH="$RESULTS"
    export NEKO_DIR="${NEKO_DIR:-$ROOT/external/neko}"
    export example="$name"
    source functions.sh
    if ! run "$name"; then
      # On NFS, rsync can encounter a transient busy .nfs file after the
      # solver has completed. The results have already been archived, so do
      # not let that cleanup-only warning abort the remaining cases.
      if grep -qE 'sender failed to remove .*\.nfs.*Device or resource busy' error.log \
          && ! grep -qE 'ERROR: (An error occurred during execution|An error occured during preparation)' error.log; then
        printf '[gauntlet] ignoring NFS cleanup warning for %s\n' "$name" >&2
        : > error.log
      else
        exit 1
      fi
    fi
  ) > "$dir/output.log" 2> "$dir/error.log"
}

run_case() {
  local name=$1
  local source_dir=$2
  local source_case=$3
  local binary=$4
  local mesh=${5:-2d}

  stage_case "$name" "$source_dir" "$source_case" "$binary" "$mesh"
  run_staged_case "$name" "$source_case"
}

copy_field() {
  local source=$1
  local destination=$2
  local field_name=$3

  mkdir -p "$destination"
  [[ -f "$source/$field_name.f00000" ]] || {
    echo "Missing field: $source/$field_name.f00000" >&2
    exit 1
  }
  cp "$source/$field_name.f00000" "$destination/$field_name.f00000"
  if [[ -f "$source/$field_name.nek5000" ]]; then
    cp "$source/$field_name.nek5000" "$destination/$field_name.nek5000"
  fi
}

run_3d_floquet_case() {
  local name=$1
  local source_dir=$2
  local source_case=$3
  local binary=$4
  local baseflow_dir=$5

  if case_has_results "$name"; then
    echo "[gauntlet] skipping $name (results already exist)"
    return 0
  fi

  stage_case "$name" "$source_dir" "$source_case" "$binary" 3d
  cp "$baseflow_dir/baseflow3d0.f00000" "$WORK/$name/"
  cp "$baseflow_dir/baseflow3d0.nek5000" "$WORK/$name/"
  cp "$SHARED/1cyl_3d.nmsh" "$WORK/$name/newton.nmsh"
  run_staged_case "$name" "$source_case"
}

prepare_2d_mesh

# Re40: steady baseflow -> direct, adjoint, and transient-growth analyses.
run_case Re40/baseflow "$CYLINDER/Re40/baseflow" 1cyl_newton.case \
  "$CYLINDER/Re40/baseflow/neko"
copy_field "$RESULTS/Re40/baseflow/viz/baseflow_newton_Re40" \
  "$WORK/Re40/baseflow" baseflow_newton_Re400

stage_case Re40/LNS "$CYLINDER/Re40/LNS" 1cyl_stability_direct.case \
  "$CYLINDER/Re40/LNS/neko"
copy_field "$RESULTS/Re40/baseflow/viz/baseflow_newton_Re40" \
  "$WORK/Re40/LNS" baseflow_newton_Re400
run_staged_case Re40/LNS 1cyl_stability_direct.case
stage_case Re40/adjoint "$CYLINDER/Re40/adjoint" 1cyl_stability_adjoint.case \
  "$CYLINDER/Re40/adjoint/neko"
copy_field "$RESULTS/Re40/baseflow/viz/baseflow_newton_Re40" \
  "$WORK/Re40/adjoint" baseflow_newton_Re400
run_staged_case Re40/adjoint 1cyl_stability_adjoint.case
stage_case Re40/transient_growth "$CYLINDER/Re40/transient_growth" \
  1cyl_transient_growth.case "$CYLINDER/Re40/transient_growth/neko"
copy_field "$RESULTS/Re40/baseflow/viz/baseflow_newton_Re40" \
  "$WORK/Re40/transient_growth" baseflow_newton_Re400
run_staged_case Re40/transient_growth 1cyl_transient_growth.case

# Re60: stabilized steady baseflow -> direct and adjoint analyses.
run_case Re60/baseflow "$CYLINDER/Re60/baseflow" 1cyl_newton.case \
  "$CYLINDER/Re60/baseflow/neko"
copy_field "$RESULTS/Re60/baseflow/viz/baseflow_newton_Re60" \
  "$WORK/Re60/baseflow/viz" baseflow_newton_Re600
run_case Re60/LNS "$CYLINDER/Re60/LNS" 1cyl_stability_direct.case \
  "$CYLINDER/Re60/LNS/neko"
run_case Re60/adjoint "$CYLINDER/Re60/adjoint" 1cyl_stability_adjoint.case \
  "$CYLINDER/Re60/adjoint/neko"

# Re180: periodic orbit -> forward and adjoint Floquet analyses.
run_case Re180/periodic_orbit "$CYLINDER/Re180/periodic_orbit" \
  1cyl_periodic_orbit.case "$CYLINDER/Re180/periodic_orbit/neko"
stage_case Re180/LNS_floquete "$CYLINDER/Re180/LNS_floquete" \
  1cyl_floquet.case "$CYLINDER/Re180/LNS_floquete/neko"
copy_field "$RESULTS/Re180/periodic_orbit/baseflow" \
  "$WORK/Re180/LNS_floquete" baseflow0
run_staged_case Re180/LNS_floquete 1cyl_floquet.case
stage_case Re180/adjoint_floquete "$CYLINDER/Re180/adjoint_floquete" \
  1cyl_adjoint_floquet.case "$CYLINDER/Re180/adjoint_floquete/neko"
copy_field "$RESULTS/Re180/periodic_orbit/baseflow" \
  "$WORK/Re180/adjoint_floquete" baseflow0
run_staged_case Re180/adjoint_floquete 1cyl_adjoint_floquet.case

# Re200: periodic orbit followed by the unstable 3D Floquet analyses.
run_case Re200/periodic_orbit "$CYLINDER/Re200/periodic_orbit" \
  1cyl_periodic_orbit.case "$CYLINDER/Re200/periodic_orbit/neko"

if [[ "$RUN_3D" == true ]]; then
  three_d="$WORK/Re180/3D"
  mkdir -p "$three_d"
  cp "$CYLINDER/Re180/LNS_floquete_3D/newton.geo" "$three_d/newton.geo"
  # Embed the requested parameters in the generated standalone prepare script.
  cat > "$three_d/prepare.sh" <<EOF
#!/usr/bin/env bash
set -euo pipefail
# Gmsh renames an existing explicit Save target instead of reliably replacing
# it.  Remove only generated 3D mesh artifacts so an interrupted preparation
# cannot make gmsh2nek consume a stale, incomplete multi-gigabyte file.
rm -f -- newton.msh newton.msh_* newton.re2 1cyl_3d.nmsh
gmsh -0 newton.geo \
  -setnumber "Spanwise domain length" "$SPAN_LENGTH" \
  -setnumber "Spanwise element layers" "$SPAN_LAYERS" >/dev/null
gmsh2nek <<'MESH_EOF' >/dev/null
3
newton
1
6 7
0 0 $SPAN_LENGTH
MESH_EOF
rea2nbin newton.re2 1cyl_3d.nmsh
EOF
  chmod u+x "$three_d/prepare.sh"
  cp "$ROOT/scripts/functions.sh" "$three_d/functions.sh"
  (
    cd "$three_d"
    export MAIN_DIR="$ROOT"
    export NEKO_DIR="${NEKO_DIR:-$ROOT/external/neko}"
    source functions.sh
    prepare
  )
  mv "$three_d/1cyl_3d.nmsh" "$SHARED/1cyl_3d.nmsh"
  python3 "$CYLINDER/Re180/LNS_floquete_3D/extrude_baseflow.py" \
    "$RESULTS/Re180/periodic_orbit/baseflow/baseflow0.f00000" \
    "$three_d/baseflow3d0.f00000" --length "$SPAN_LENGTH" --layers "$SPAN_LAYERS"
  printf '%s %s\n' "$SPAN_LENGTH" "$SPAN_LAYERS" > "$SPAN_CONFIG_FILE"

  run_3d_floquet_case Re180/LNS_floquete_3D \
    "$CYLINDER/Re180/LNS_floquete_3D" 1cyl_floquet_3D.case \
    "$CYLINDER/Re180/LNS_floquete_3D/neko" "$three_d"
  run_3d_floquet_case Re180/adjoint_floquete_3D \
    "$CYLINDER/Re180/adjoint_floquete_3D" 1cyl_adjoint_floquet_3D.case \
    "$CYLINDER/Re180/adjoint_floquete_3D/neko" "$three_d"

  # Reuse the same spanwise mesh for Re200, but extrude its independently
  # converged periodic orbit and stage both 3D Floquet operators at Re=200.
  three_d_200="$WORK/Re200/3D"
  mkdir -p "$three_d_200"
  python3 "$CYLINDER/Re200/LNS_floquete_3D/extrude_baseflow.py" \
    "$RESULTS/Re200/periodic_orbit/baseflow/baseflow0.f00000" \
    "$three_d_200/baseflow3d0.f00000" \
    --length "$SPAN_LENGTH" --layers "$SPAN_LAYERS"

  run_3d_floquet_case Re200/LNS_floquete_3D \
    "$CYLINDER/Re200/LNS_floquete_3D" 1cyl_floquet_3D.case \
    "$CYLINDER/Re200/LNS_floquete_3D/neko" "$three_d_200"
  run_3d_floquet_case Re200/adjoint_floquete_3D \
    "$CYLINDER/Re200/adjoint_floquete_3D" 1cyl_adjoint_floquet_3D.case \
    "$CYLINDER/Re200/adjoint_floquete_3D/neko" "$three_d_200"
fi

echo "[gauntlet] complete: $GAUNTLET"
