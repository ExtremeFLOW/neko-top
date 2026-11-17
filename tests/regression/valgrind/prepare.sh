#!/usr/bin/bash

# Handle options
Nx=16 && Ny=8 && Nz=8

# ============================================================================ #
# Ensure Neko can be found and set default mesh size

if [ -n "$NEKO_DIR" ]; then
    PATH=$NEKO_DIR/bin:$PATH
fi

if command -v genmeshbox >/dev/null 2>&1; then
    echo -e "Neko tool 'genmeshbox' not found." >&2
    echo -e "Please ensure Neko is installed and in your PATH." >&2
    echo -e "Alternatively, set the NEKO_DIR environment variable." >&2
    exit 1
fi

# ============================================================================ #
# Generate mesh and run case

echo "Generating mesh with dimensions: $Nx $Ny $Nz"
if [ ! -f "box.nmsh" ]; then
    genmeshbox 0 4 0 1 0 1 $Nx $Ny $Nz .false. .false. .false.
else
    echo "Mesh file 'box.nmsh' already exists. Skipping mesh generation."
fi

# End of file
# ============================================================================ #
