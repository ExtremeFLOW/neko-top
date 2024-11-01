Nx=20 && Ny=20 && Nz=1
echo "Generating mesh with dimensions: $Nx $Ny $Nz"
genmeshbox 0 1 0 1 0 1 $Nx $Ny $Nz .false. .true. .true.
