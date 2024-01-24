import meshio

# Read the mesh
mesh = meshio.read("data_local/output_00059.vtu")
meshio.write("data_local/output_00059.vtk", mesh)
