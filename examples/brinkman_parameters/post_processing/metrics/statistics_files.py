# Imports for Bertie tools
import os
os.environ["PYNEKTOOLS_HIDE_LOG"] = 'true'
from mpi4py import MPI #equivalent to the use of MPI_init() in C
from pynektools.io.ppymech.neksuite import pynekread
from pynektools.datatypes.msh import Mesh
from pynektools.datatypes.field import FieldRegistry
import numpy as np
import pickle
import matplotlib.pyplot as plt

def wake(file_name: str, plot_params, cache_dir: str = None) -> dict:

    if cache_dir is not None:
        file_ext = os.path.splitext(file_name)[1]
        cache_file = os.path.join(
            cache_dir, "wake_lines",
            os.path.relpath(file_name,
                            cache_dir).replace("../", "").replace("/", "_"))

        cache_file = cache_file + ".pkl"
        cache_exists = os.path.exists(cache_file)
        file_exists = os.path.exists(file_name)

        # If the cache file exists and is not older than the file, load it
        mod_time = os.path.getmtime(file_name) if file_exists else 0
        cache_time = os.path.getmtime(cache_file) if cache_exists else 0
        if cache_time > mod_time:
            with open(cache_file, "rb") as f:
                return pickle.load(f)

    # ------------------------------------------------------------------------ #
    # Main computation
    name = file_name + '/fluid_stats0/fluid_stats0'
    case = {}
    mesh, fld = load_stats_file(name, 2)
    y_lims = plot_params["wake_y_lim"]
    n_pts = plot_params["wake_n_pts"]
    case["wake_positions"] = plot_params["wake_positions"]
    case["wake_xyz"], case["wake_u"]  = generate_wake_lines(mesh, fld, case["wake_positions"], y_lims, n_pts)


    # ------------------------------------------------------------------------ #
    # Caching the results

    if cache_dir is not None:
        os.makedirs(os.path.dirname(cache_file), exist_ok=True)
        with open(cache_file, "wb") as f:
            pickle.dump(case, f)

    return case

def load_stats_file(file_name, fld_number):

    # Import required modules
    from mpi4py import MPI #equivalent to the use of MPI_init() in C
    
    # Get mpi info
    comm = MPI.COMM_WORLD
    
    # Hide the log for the notebook. Not recommended when running in clusters as it is better you see what happens
    #import os
    #os.environ["PYNEKTOOLS_HIDE_LOG"] = 'true'
    
    from pynektools.io.ppymech.neksuite import pynekread
    from pynektools.datatypes.msh import Mesh
    

    # set mesh and fld
    
    msh = Mesh(comm, create_connectivity=True)
    fld = FieldRegistry(comm)
    
    # remember the mesh is only stored in the first one!
    # here we're just getting the mesh
    pynekread(file_name+'.f00000', comm, data_dtype=np.double, msh=msh)
    
    # TODO
    # Right now I'm just playing around with one, but ultimately we will take an "average of the averages"
    pynekread(file_name + '.f' + "{:05d}".format(fld_number), comm, data_dtype=np.double, fld=fld)

    return msh, fld

def generate_wake_lines(msh, fld, wake_positions, y_lims, n_pts):
    # Import required modules
    from mpi4py import MPI #equivalent to the use of MPI_init() in C
    import matplotlib.pyplot as plt
    import numpy as np

    # Get mpi info
    comm = MPI.COMM_WORLD

    # Hide the log for the notebook. Not recommended when running in clusters as it is better you see what happens
    import os
    os.environ["PYNEKTOOLS_HIDE_LOG"] = 'true'

    from pynektools.io.ppymech.neksuite import pynekread
    from pynektools.datatypes.msh import Mesh
    from pynektools.datatypes.coef import Coef
    from pynektools.datatypes.field import FieldRegistry
    from pynektools.postprocessing.file_indexing import index_files_from_folder

    from pynektools.interpolation.probes import Probes
    import pynektools.interpolation.utils as interp_utils
    import pynektools.interpolation.pointclouds as pcs

    # set up the probe locations
    if comm.Get_rank() == 0 :
        # Wake line
        y_1d = pcs.generate_1d_arrays(y_lims, n_pts, mode="equal")
        x_1d = wake_positions
        z_1d = [0.5]
        x, y, z = np.meshgrid(x_1d, y_1d, z_1d, indexing='ij')
        # Write the points for future use
        xyz = interp_utils.transform_from_array_to_list(len(wake_positions),n_pts,1,[x, y, z])
    else:
        xyz = None


    # This sets up the probes interpolation
    probes = Probes(comm, probes = xyz, msh = msh, point_interpolator_type='multiple_point_legendre_numpy', \
                             max_pts=256, find_points_comm_pattern='point_to_point')
    # here we actually probe
    probes.interpolate_from_field_list(0, [fld.registry['u']], comm, write_data=False)
    # package them nicely
    int_fields = interp_utils.transform_from_list_to_array(len(wake_positions),n_pts,1,probes.interpolated_fields)
    wake_u = np.copy(int_fields[1])
    # I think we only really need u...
    return xyz, wake_u

def init_plot_a_wake_line(plot_params):
    # gotta assume they're all the same
    wake_positions = plot_params["wake_positions"]
    fig, ax = plt.subplots(1, len(wake_positions), figsize=(10, 5), dpi = 200)
    return fig, ax

def plot_a_wake_line(case, plot_params, ax, fig, color, linestyle):
    # gotta assume they're all the same
    wake_positions = plot_params["wake_positions"]
    for i in range(len(wake_positions)):
        wake_index = np.where(case["wake_xyz"][:,0] == wake_positions[i])
        ax[i].plot(case["wake_u"][i,:].flatten(), case["wake_xyz"][wake_index, 1].flatten(), label=case["name"])

def finalize_plot_a_wake_line(plot_params, ax, fig, output_filename):
    wake_positions = plot_params["wake_positions"]
    U_lims = plot_params["wake_U_lim"]
    y_lims = plot_params["wake_y_lim"]
    for i in range(len(wake_positions)):
        ax[i].set_title('Wake at x ='+str(wake_positions[i]))
        ax[i].set_xlim(U_lims)
        ax[i].set_ylim(y_lims)
    ax[len(wake_positions)-1].legend()
    fig.savefig(output_filename)

def calculate_forces(msh, fld, ring_radii, n_points, Re):
    from mpi4py import MPI #equivalent to the use of MPI_init() in C
    from pynektools.io.ppymech.neksuite import pynekread
    from pynektools.datatypes.msh import Mesh
    from pynektools.datatypes.coef import Coef
    from pynektools.datatypes.field import FieldRegistry

    # Hide the log for the notebook. Not recommended when running in clusters as it is better you see what happens
    import os
    os.environ["PYNEKTOOLS_HIDE_LOG"] = 'true'
    
    from pynektools.io.ppymech.neksuite import pynekread
    from pynektools.datatypes.msh import Mesh
    from pynektools.datatypes.coef import Coef
    from pynektools.datatypes.field import FieldRegistry
    from pynektools.postprocessing.file_indexing import index_files_from_folder

    from pynektools.interpolation.probes import Probes
    import pynektools.interpolation.utils as interp_utils
    import pynektools.interpolation.pointclouds as pcs

    # Get mpi info
    comm = MPI.COMM_WORLD
    
    # compute derivatives
    coef = Coef(msh, comm, get_area=False)
    dudx = coef.dudxyz(fld.registry['u'], coef.drdx, coef.dsdx, coef.dtdx)
    dudy = coef.dudxyz(fld.registry['u'], coef.drdy, coef.dsdy, coef.dtdx)
    dvdx = coef.dudxyz(fld.registry['v'], coef.drdx, coef.dsdx, coef.dtdx)
    dvdy = coef.dudxyz(fld.registry['v'], coef.drdy, coef.dsdy, coef.dtdx)

    # set up probes
    if comm.Get_rank() == 0 :
    
        # Generate the 1D mesh
        r_1d = np.array(ring_radii)
        th_1d = pcs.generate_1d_arrays([0, 2*np.pi], n_points, mode="equal")
        z_1d = np.array([0.5])
    
        # Generate a 3D mesh
        r, th, z = np.meshgrid(r_1d, th_1d, z_1d, indexing='ij')
        x = r*np.cos(th)
        y = r*np.sin(th)
    
        # Array the points as a list of probes
        xyz = interp_utils.transform_from_array_to_list(len(ring_radii),n_points,1,[x, y, z])
    
    else:
        xyz = None
    
    
    # This sets up the probes interpolation
    ring_probes = Probes(comm, probes = xyz, msh = msh, point_interpolator_type='multiple_point_legendre_numpy', \
                             max_pts=256, find_points_comm_pattern='point_to_point', output_fname='/tmp/probes.csv')
    
    # here we actually probe
    ring_probes.interpolate_from_field_list(0, [dudx], comm, write_data=False)
    int_fields = interp_utils.transform_from_list_to_array(len(ring_radii),n_points,1,ring_probes.interpolated_fields)
    dudx = np.copy(int_fields[1])
    
    ring_probes.interpolate_from_field_list(0, [dudy], comm, write_data=False)
    int_fields = interp_utils.transform_from_list_to_array(len(ring_radii),n_points,1,ring_probes.interpolated_fields)
    dudy = np.copy(int_fields[1])
    
    ring_probes.interpolate_from_field_list(0, [dvdx], comm, write_data=False)
    int_fields = interp_utils.transform_from_list_to_array(len(ring_radii),n_points,1,ring_probes.interpolated_fields)
    dvdx = np.copy(int_fields[1])
    
    ring_probes.interpolate_from_field_list(0, [dvdy], comm, write_data=False)
    int_fields = interp_utils.transform_from_list_to_array(len(ring_radii),n_points,1,ring_probes.interpolated_fields)
    dvdy = np.copy(int_fields[1])
    
    ring_probes.interpolate_from_field_list(0, [fld.registry['p']], comm, write_data=False)
    int_fields = interp_utils.transform_from_list_to_array(len(ring_radii),n_points,1,ring_probes.interpolated_fields)
    p = np.copy(int_fields[1])
    
    f_x_visc = np.zeros((len(ring_radii), n_points))
    f_y_visc = np.zeros((len(ring_radii), n_points))
    f_x_p    = np.zeros((len(ring_radii), n_points))
    f_y_p    = np.zeros((len(ring_radii), n_points))
    f_x_tot  = np.zeros((len(ring_radii), n_points))
    f_y_tot  = np.zeros((len(ring_radii), n_points))
    
    f_x_visc_int = np.zeros((len(ring_radii)))
    f_y_visc_int = np.zeros((len(ring_radii)))
    f_x_p_int    = np.zeros((len(ring_radii)))
    f_y_p_int    = np.zeros((len(ring_radii)))
    f_x_tot_int  = np.zeros((len(ring_radii)))
    f_y_tot_int  = np.zeros((len(ring_radii)))
    
    for i in range(len(ring_radii)):
        for j in range(n_points):
            # viscous
            f_x_visc[i,j] = ( dudx[i,j] * 2          * np.cos(th_1d[j]) + \
                             (dudy[i,j] + dvdx[i,j]) * np.sin(th_1d[j]))/Re
            f_y_visc[i,j] = ((dudy[i,j] + dvdx[i,j]) * np.cos(th_1d[j]) + \
                              dvdy[i,j] * 2          * np.sin(th_1d[j]))/Re
            
            #pressure
            f_x_p[i,j] = -p[i,j] * np.cos(th_1d[j])
            f_y_p[i,j] = -p[i,j] * np.sin(th_1d[j])

            # total
            f_x_tot[i,j] = f_x_visc[i,j] + f_x_p[i,j]
            f_y_tot[i,j] = f_y_visc[i,j] + f_y_p[i,j]

    # Here we integrate them around to get a single number
    for i in range(len(ring_radii)):
        ri1 = (n_points)*i
        ri2 = (n_points )*(i+1)
        dS = (2 * np.pi * ring_radii[i]) / n_points
        
        f_x_visc_int[i] = np.sum(f_x_visc[i,:]) * dS
        f_y_visc_int[i] = np.sum(f_y_visc[i,:]) * dS
        f_x_p_int[i]    = np.sum(f_x_p[i,:]) * dS
        f_y_p_int[i]    = np.sum(f_y_p[i,:]) * dS
        
        f_x_tot_int[i] = f_x_visc_int[i] + f_x_p_int[i]
        f_y_tot_int[i] = f_y_visc_int[i] + f_y_p_int[i]

        force_measure = {}
        force_measure["fx_tot"] = f_x_tot
        force_measure["fx_p"] = f_x_p
        force_measure["fx_visc"] = f_x_visc
        force_measure["fy_tot"] = f_y_tot
        force_measure["fy_p"] = f_y_p
        force_measure["fy_visc"] = f_y_visc
        force_measure["fx_tot_int"] = f_x_tot_int
        force_measure["fx_p_int"] = f_x_p_int
        force_measure["fx_visc_int"] = f_x_visc_int
        force_measure["fy_tot_int"] = f_y_tot_int
        force_measure["fy_p_int"] = f_y_p_int
        force_measure["fy_visc_int"] = f_y_visc_int
        force_measure["ring"] = ring_radii

    return force_measure

def plot_ring(case_list, rings_plot, output_filename):
    import pynektools.interpolation.pointclouds as pcs
    # assume they're all the same
    n_points = case_list[0]["stats_forces"]["fx_tot"].shape[1]
    th_1d = pcs.generate_1d_arrays([0, 2*np.pi], n_points, mode="equal")

    fig, ax = plt.subplots(len(rings_plot),1, figsize=(10, 5*len(rings_plot)), dpi = 200)

    # lets just look at 2
    for i in range(len(rings_plot)):
        for j in range(len(case_list)):
            case = case_list[j]
            ring_index = np.where(np.array(case["stats_forces"]["ring"]) == np.array(rings_plot[i]))
            F_tot = np.sqrt(case["stats_forces"]["fx_tot"][ring_index,:].flatten()**2 + case["stats_forces"]["fy_tot"][ring_index,:].flatten()**2)
            ax[i].plot(th_1d/np.pi*180, F_tot, label=case["name"])
        ax[i].set_title('Ring at r ='+str(rings_plot[i]))
        ax[i].set_xlabel(r'$\theta$')
        ax[i].set_ylabel(r'$|F|$')
    ax[0].legend()
    plt.savefig(output_filename)

def plot_total_against_radius(case_list, output_file):
    # Maybe look at how lift and drag changes with the circle we take:
    fig, ax = plt.subplots(1,2, figsize=(10, 5), dpi = 200)
    for j in range(len(case_list)):
        case = case_list[j]
        ax[0].plot(case["stats_forces"]["ring"], case["stats_forces"]["fx_tot_int"], label=case["name"])
        ax[1].plot(case["stats_forces"]["ring"], case["stats_forces"]["fy_tot_int"], label=case["name"])
    ax[0].set_title('Drag')
    ax[0].set_xlabel('r')
    ax[1].set_title('Lift')
    ax[1].set_xlabel('r')
    ax[1].legend()
    plt.savefig(output_file)
