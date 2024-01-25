import meshio

# Read the mesh
mesh = meshio.read("data_local/petsc_static_mixer_lowest.vtu")
meshio.write("data_local/petsc_static_mixer_lowest.vtk", mesh, binary=False)
