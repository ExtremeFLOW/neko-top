def calc_lift_and_drag(root_name,case_name,method,Re):
    # somehow search for all the rings! If i was better at python I would use the os somehow to find everything starting with 'circ'
    # but I'm not.. so I'm listing them manually.
    circ_list = ['048', '050', '0501', '0502', '0505', '051', '052', '055']
    # First is lift and drag
    if method == 'meshed':
        # here we just need the logfile
        path_name = root_name + case_name + '/cylinder.log' 
        t, fx_tot, fx_p, fx_visc, fy_tot, fy_p, fy_visc = read_force_torque(path_name)
        force_measure = {}
        force_measure["t"] = t
        force_measure["fx_tot"] = fx_tot
        force_measure["fx_p"] = fx_p
        force_measure["fx_visc"] = fx_visc
        force_measure["fy_tot"] = fy_tot
        force_measure["fy_p"] = fy_p
        force_measure["fy_visc"] = fy_visc
        force_measure["type"] = "meshed"
        all_forces = [force_measure]
        # We could in principal use method 2, but it's not required I guess.
    elif method == 'brinkman':
        # Now we have concentric rings
        all_forces = []
        for i in range(len(circ_list)):
            path_name = root_name + case_name + '/circ_' + circ_list[i] + '.csv'
            t, fx_tot, fx_p, fx_visc, fy_tot, fy_p, fy_visc = surface_integral_from_probes(path_name, Re)
            force_measure = {}
            force_measure["t"] = t
            force_measure["fx_tot"] = fx_tot
            force_measure["fx_p"] = fx_p
            force_measure["fx_visc"] = fx_visc
            force_measure["fy_tot"] = fy_tot
            force_measure["fy_p"] = fy_p
            force_measure["fy_visc"] = fy_visc
            force_measure["type"] = 'circ_' + circ_list[i] + '.csv'
            all_forces.append(force_measure)

        # we also have the 3rd measure
        path_name = root_name + case_name + "/cylinder.log"
        t, fx_tot, fy_tot = read_brinkman_force(path_name)
        # let's just fill the others with zeros to indicate this is the strange one:
        fx_p = np.zeros(len(t))
        fx_visc = np.zeros(len(t))
        fy_p = np.zeros(len(t))
        fy_visc = np.zeros(len(t))
        force_measure = {}
        force_measure["t"] = t
        force_measure["fx_tot"] = fx_tot
        force_measure["fx_p"] = fx_p
        force_measure["fx_visc"] = fx_visc
        force_measure["fy_tot"] = fy_tot
        force_measure["fy_p"] = fy_p
        force_measure["fy_visc"] = fy_visc
        force_measure["type"] = 'method3'
        all_forces.append(force_measure)
    elif method == 'idw':
        # I guess we do the same as brinkman without the last step
        # Now we have concentric rings
        all_forces = []
        for i in range(len(circ_list)):
            path_name = root_name + case_name + '/circ_' + circ_list[i] + '.csv'
            t, fx_tot, fx_p, fx_visc, fy_tot, fy_p, fy_visc = surface_integral_from_probes(path_name, Re)
            force_measure = {}
            force_measure["t"] = t
            force_measure["fx_tot"] = fx_tot
            force_measure["fx_p"] = fx_p
            force_measure["fx_visc"] = fx_visc
            force_measure["fy_tot"] = fy_tot
            force_measure["fy_p"] = fy_p
            force_measure["fy_visc"] = fy_visc
            force_measure["type"] = 'circ_' + circ_list[i] + '.csv'
            all_forces.append(force_measure)
    else:
        print('INCORRECT METHOD!')

    return all_forces

def plot_lift_drag(case_list, lift_axis, drag_axis):
    fig, ax = plt.subplots(2,1, figsize=(10, 5), dpi = 200)
    for j in range(len(case_list)):
        case = case_list[j]
        for i in range(len(case["forces"])):
            force_measure = case["forces"][i]
            # I can tell already we only only want 0.5
            # but this can be changed
            if force_measure["type"] == 'circ_050' or force_measure["type"] == 'meshed' or force_measure["type"] == 'method3':
                ax[0].plot(force_measure["t"], force_measure["fy_tot"],label=case["name"] + force_measure["type"])
                ax[1].plot(force_measure["t"], force_measure["fx_tot"], label=case["name"] + force_measure["type"])
    ax[0].set_ylabel("Lift")
    ax[1].set_ylabel("Drag")
    ax[1].set_xlabel("Time")
    ax[0].set_ylim(lift_axis)
    ax[1].set_ylim(drag_axis)
    ax[1].legend()
    plt.show()
    
def read_force_torque(path_name):
    # note!
    # I bet this isn't the best for long logfiles!
    t_list = []
    fx_tot_list = []
    fx_visc_list = []
    fx_p_list = []
    fy_tot_list = []
    fy_visc_list = []
    fy_p_list = []
    
    with open(path_name) as file:
        for line in file:
            s = line.rstrip()
            index = s.find("forcex")
            if index >= 1: 
                numbers = np.array(line.rstrip().split())
                t_list.append(numbers[1].astype(float))
                fx_tot_list.append(numbers[2].astype(float))
                fx_p_list.append(numbers[3].astype(float))
                fx_visc_list.append(float(numbers[4][:-1]))
            index = s.find("forcey")
            if index >= 1:
                numbers = np.array(line.rstrip().split())
                fy_tot_list.append(numbers[2].astype(float))
                fy_p_list.append(numbers[3].astype(float))
                fy_visc_list.append(float(numbers[4][:-1]))
            
    t = np.array(t_list)
    fx_tot = np.array(fx_tot_list)
    fx_p = np.array(fx_p_list)
    fx_visc = np.array(fx_visc_list)
    fy_tot = np.array(fy_tot_list)
    fy_p = np.array(fy_p_list)
    fy_visc = np.array(fy_visc_list)

    return t, fx_tot, fx_p, fx_visc, fy_tot, fy_p, fy_visc



def surface_integral_from_probes(input_filename, Re):
    # Read the first line to get the number of points
    with open(input_filename, 'r') as file:
        first_line = file.readline().strip()
        n = int(first_line.split(',')[0])  # Extract the first value as the number of points

    # Read coordinates
    coordinates = pd.read_csv(input_filename, skiprows=1, nrows=n, header=None)
    x_coords = coordinates.values[:, 0]  # All x-coordinates (first column)
    y_coords = coordinates.values[:, 1]  # All y-coordinates (second column)

    # Read the data from row n+1 onwards
    data = pd.read_csv(input_filename, skiprows=n + 1, header=None)



    # Extract time and variable columns
    time = data[0]  # Column A is time
    u = data[1] # Column B is u
    v = data[2] # Column C is v
    w = data[3] # Column D is w
    p = data[4] # Column E is p
    # print ("data.shape[1]=" , data.shape[1])
    if data.shape[1] > 8:
        dudx = data[5]     # Column F is dudx
        dudy = data[6]     # Column G is dudy
        dvdx = data[7]     # Column H is dvdx
        dvdy = data[8]     # Column I is dvdy
        derivative_available = True
    print("Loading data from file (",os.path.basename(input_filename),") successfully done. n= ", n, " points/probs detected.")
    print(derivative_available)
    radius = np.sqrt(x_coords[0]**2 + y_coords[0]**2)  # Radius of the circle
    # Calculate theta for each probe based on x and y coordinates
    theta = np.arctan2(y_coords, x_coords)  # Angle in radians for each probe



    num_timesteps = int(len(time) / n)
    time_reshaped = np.zeros((n, num_timesteps))
    p_reshaped = np.zeros((n, num_timesteps))
    u_reshaped = np.zeros((n, num_timesteps))
    v_reshaped = np.zeros((n, num_timesteps))

    dudx_reshaped = np.zeros((n, num_timesteps))
    dudy_reshaped = np.zeros((n, num_timesteps))
    dvdx_reshaped = np.zeros((n, num_timesteps))
    dvdy_reshaped = np.zeros((n, num_timesteps))

    # Loop to manually fill the reshaped arrays
    for i in range(0, n-1, 1):
        # For each probe i, fill the corresponding row in each reshaped array
        time_reshaped[i,:] = time.values[i::n]
        p_reshaped[i,:] = p.values[i::n]
        u_reshaped[i,:]= u.values[i::n]
        v_reshaped[i,:] = v.values[i::n]

        dudx_reshaped[i,:] = dudx.values[i::n]
        dudy_reshaped[i,:] = dudy.values[i::n]
        dvdx_reshaped[i,:] = dvdx.values[i::n]
        dvdy_reshaped[i,:] = dvdy.values[i::n]

    # Initialize arrays to store drag, lift, and shear components
    drag = np.zeros(num_timesteps)
    lift = np.zeros(num_timesteps)
    shear_drag = np.zeros(num_timesteps)
    shear_lift = np.zeros(num_timesteps)


    # Calculate the arc length per segment

    dS = (2 * np.pi * radius) / n


    # Loop over each time step to calculate drag and lift
    for t in range(num_timesteps):
        # Get pressure and angle for each probe at this time step
        p_t = p_reshaped[:, t]



        # Get shear components (du/dn and dv/dn) for each probe at this time step
        #du_dn = np.cos(theta) * dudx_reshaped[:, t] + np.sin(theta) * dudy_reshaped[:, t]
        #dv_dn = np.cos(theta) * dvdx_reshaped[:, t] + np.sin(theta) * dvdy_reshaped[:, t]

        # Compute the shear drag and lift components (using du/dn and dv/dn)
        #shear_drag_components = du_dn * np.cos(theta)  # Shear contribution to drag
        #shear_lift_components = dv_dn * np.sin(theta)  # Shear contribution to lift

        shear_drag_components = ( dudx_reshaped[:, t] * 2                    * np.cos(theta) + \
                                 (dvdx_reshaped[:, t] + dudy_reshaped[:, t]) * np.sin(theta))/Re  # Shear contribution to drag
        shear_lift_components = ((dvdx_reshaped[:, t] + dudy_reshaped[:, t]) * np.cos(theta) + \
                                  dvdy_reshaped[:, t] * 2                    * np.sin(theta))/Re  # Shear contribution to lift

        # Compute the drag and lift components from pressure
        drag_components = -p_t * np.cos(theta)
        lift_components = -p_t * np.sin(theta)

        # Integrate around the circle (sum and multiply by 2 * pi / n for full circle)
        drag[t] = np.sum(drag_components) * dS
        lift[t] = np.sum(lift_components) * dS
        shear_drag[t] = np.sum(shear_drag_components) * dS
        shear_lift[t] = np.sum(shear_lift_components) * dS

    # Total drag and lift (pressure + shear)(1/Re)(Re=200)
    total_drag = drag + shear_drag
    total_lift = lift + shear_lift
    time_values = time_reshaped[0,:] # Extract the time from the last probe's time array

    return time_values, total_drag, drag, shear_drag, total_lift, lift, shear_lift


def read_brinkman_force(path_name):
    t_list = []
    fx_list = []
    fy_list = []

    with open(path_name) as file:
        for line in file:
            s = line.rstrip()
            index = s.find("Lift")
            if index >= 1:
                numbers = np.array(line.rstrip().split())
                t_list.append(numbers[4].astype(float))
                fy_list.append(numbers[2].astype(float))
            index = s.find("Drag")
            if index >= 1:
                numbers = np.array(line.rstrip().split())
                fx_list.append(numbers[2].astype(float))

    t = np.array(t_list)
    fx = np.array(fx_list)
    fy = np.array(fy_list)
    return t, fx, fy

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

def plot_wake_lines(case_list):
    # gotta assume they're all the same
    wake_positions = case_list[0]["wake_positions"]
    fig, ax = plt.subplots(1, len(wake_positions), figsize=(10, 5), dpi = 200)
    for i in range(len(wake_positions)):
        for j in range(len(case_list)):
            case = case_list[j]
            wake_index = np.where(case["wake_xyz"][:,0] == wake_positions[i])
            ax[i].plot(case["wake_u"][i,:].flatten(), case["wake_xyz"][wake_index, 1].flatten(), label=case["name"])
        ax[i].set_title('Wake at x ='+str(wake_positions[i]))
        ax[i].set_xlim(U_lims)
        ax[i].set_ylim(y_lims)
    ax[len(wake_positions)-1].legend()
    plt.show()

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
                             max_pts=256, find_points_comm_pattern='point_to_point')
    
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

def plot_ring(case_list, rings_plot):
    th_1d = pcs.generate_1d_arrays([0, 2*np.pi], n_points, mode="equal")

    fig, ax = plt.subplots(len(rings_plot),1, figsize=(10, 5*len(rings_plot)), dpi = 200)

    # lets just look at 2
    for i in range(len(rings_plot)):
        for j in range(len(case_list)):
            case = case_list[j]
            ring_index = np.where(np.array(case["stats_forces"]["ring"]) == np.array(rings_plot[i]))
            print(ring_index)
            F_tot = np.sqrt(case["stats_forces"]["fx_tot"][ring_index,:].flatten()**2 + case["stats_forces"]["fy_tot"][ring_index,:].flatten()**2)
            ax[i].plot(th_1d/np.pi*180, F_tot, label=case["name"])
        ax[i].set_title('Ring at r ='+str(np.array(ring_radii)[ring_index]))
        ax[i].set_xlabel(r'$\theta$')
        ax[i].set_ylabel(r'$|F|$')
    ax[0].legend()
    plt.show()

def plot_total_against_radius(case_list):
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
    plt.show()
