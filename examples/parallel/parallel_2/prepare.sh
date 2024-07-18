#!/usr/bin/bash

# ============================================================================ #
# Ensure Neko can be found and set default mesh size

if [ "$NEKO_DIR" ]; then
    PATH=$NEKO_DIR/bin:$PATH
fi

if [[ -z $(which genmeshbox) ]]; then
    echo -e "Neko not found." >&2
    echo -e "Please ensure Neko is installed and in your PATH." >&2
    echo -e "Alternatively, set the NEKO_DIR environment variable." >&2
    exit 1
fi

# ============================================================================ #
# Generate mesh and run case
Nx=32 && Ny=8 && Nz=8

echo "Generating mesh with dimensions: $Nx $Ny $Nz"
genmeshbox 0 4 0 1 0 1 $Nx $Ny $Nz .false. .false. .false.
prepart box.nmsh 2

# End of file
# ============================================================================ #
