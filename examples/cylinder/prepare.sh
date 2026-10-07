#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo=$(cd "$root/../.." && pwd)

lib_path="$repo/external/json-fortran/lib"
if [[ -n "${HDF5_ROOT:-}" ]]; then
  lib_path="$lib_path:$HDF5_ROOT/lib:$HDF5_ROOT/lib64"
fi
if [[ -d "$repo/external/hdf5/lib" ]]; then
  lib_path="$lib_path:$repo/external/hdf5/lib"
fi
if [[ -d "$repo/external/hdf5/lib64" ]]; then
  lib_path="$lib_path:$repo/external/hdf5/lib64"
fi
export LD_LIBRARY_PATH="$lib_path:${LD_LIBRARY_PATH:-}"

rea2nbin_bin=${REA2NBIN:-"$repo/external/neko/bin/rea2nbin"}
"$rea2nbin_bin" "$root/_shared/1cyl.re2" "$root/_shared/1cyl.nmsh"

echo "Wrote $root/_shared/1cyl.nmsh"
