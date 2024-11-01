Nx=80 && Ny=40 && Nz=1
echo "Generating mesh with dimensions: $Nx $Ny $Nz"
genmeshbox -15 35 -15 15 0 1 $Nx $Ny $Nz .false. .false. .true.
