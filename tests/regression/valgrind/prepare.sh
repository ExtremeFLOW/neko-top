#!/usr/bin/bash

# Assign test size
Nx=${VALGRIND_MESH_NX:-16}
Ny=${VALGRIND_MESH_NY:-8}
Nz=${VALGRIND_MESH_NZ:-8}

# ============================================================================ #
# Ensure Neko can be found and set default mesh size

if [ -n "$NEKO_DIR" ]; then
    PATH=$NEKO_DIR/bin:$PATH
fi

if ! command -v genmeshbox >/dev/null 2>&1; then
    echo >&2 "Neko tool 'genmeshbox' not found."
    echo >&2 "Please ensure Neko is installed and in your PATH."
    echo >&2 "Alternatively, set the NEKO_DIR environment variable."
    exit 1
fi

# ============================================================================ #
# Generate mesh and run case


echo "Generating mesh with dimensions: $Nx $Ny $Nz"
genmeshbox 0 4 0 1 0 1 $Nx $Ny $Nz .false. .false. .false.


# End of file
# ============================================================================ #
