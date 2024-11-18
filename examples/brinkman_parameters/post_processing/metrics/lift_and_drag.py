""""Compute the lift and drag of the flow"""

import os
import pickle
import numpy as np
from pynektools.io.read_probes import ProbesReader

# ============================================================================ #
# Main function for computing the separation angle


def surface_integral_lift_and_drag(file_name: str, Re: float, cache_dir: str = None) -> dict:
    """
    For a given file with a given Re number, computes lift and drag for the probes defined in the file

    Parameters
    ----------
    file_name : str
        The name of the file to run the benchmark on.
    Re : float
        Re number for this case. (can also be extracted from the file_name)
    cache_dir : str, optional
        The directory to store the cache files. Cached files will be used if
        they exist. If None, no caching will be used. Default is None.

    Returns
    -------
    results : dict
        The results of the benchmark.
        t: np.array
            The time steps.
        fx_tot: np.ndarray
            Total_drag vs time
        fx_p: np.array
            Pressure drag vs time.
        fx_visc: np.array
            Shear drag vs time.
        fy_tot: np.ndarray
            Total_lift vs time
        fy_p: np.array
            Pressure lift vs time.
        fy_visc: np.array
            Shear_lift vs time.
    """
    # ------------------------------------------------------------------------ #
    # Caching
    #
    # We are going to cache the results of the benchmark to avoid recomputing
    # the results. This is useful when running the benchmark multiple times. The
    # cache file will be stored in the cache_dir directory. The benchmark
    # results will be cached using pickle and be named after the file_name
    # relative to the cache_dir.

    if cache_dir is not None:
        file_ext = os.path.splitext(file_name)[1]
        cache_file = os.path.join(
            cache_dir, "lift_and_drag",
            os.path.relpath(file_name,
                            cache_dir).replace("../", "").replace("/", "_"))
        cache_file = cache_file.replace(file_ext, ".pkl")

        if os.path.exists(cache_file):
            with open(cache_file, "rb") as f:
                return pickle.load(f)

    # ------------------------------------------------------------------------ #
    # Main computation

    # Read in the file and setup the data
    probes = ProbesReader(file_name)

    points = probes.points
    times = probes.times

    fields = np.asarray((probes.fields["p"],probes.fields["u"], probes.fields["v"],probes.fields["du_dx"], probes.fields["du_dy"],probes.fields["dv_dx"], probes.fields["dv_dy"]))
    del probes
    
    # just use the coordinates of the first prob to compute the radius for all probs
    # note that we consider all probs are on a circle
    radius = np.sqrt(points[0,0]**2 + points[0,1]**2)  # Radius of the circle
    theta = np.arctan2(points[:,1], points[:,0])  # Angle in radians for each probe
    num_timesteps = times.size

    # print(points.shape)
    # print(times.shape)
    # print(fields[0])
    # print(times.size)
    # Initialize arrays to store drag, lift, and shear components
    drag = np.zeros(num_timesteps)
    lift = np.zeros(num_timesteps)
    shear_drag = np.zeros(num_timesteps)
    shear_lift = np.zeros(num_timesteps)
    
    # Calculate the arc length per segment
    dS = (2 * np.pi * radius) / radius.size
    
    

    # Loop over the probs for each time step to calculate drag and lift for that time
    for t in range(num_timesteps):
        # Get pressure and angle for each probe at this time step
        # note that p = fields[0], u = fields[1], v = fields[2], du_dx = fields[3]
        # du_dy = fields[4], dv_dx = fields[5], dv_dy = fields[6]

        '''
        shear_drag_components = ( dudx_reshaped[:, t] * 2                    * np.cos(theta) + \
                                 (dvdx_reshaped[:, t] + dudy_reshaped[:, t]) * np.sin(theta))/Re  # Shear contribution to drag
        shear_lift_components = ((dvdx_reshaped[:, t] + dudy_reshaped[:, t]) * np.cos(theta) + \
                                  dvdy_reshaped[:, t] * 2                    * np.sin(theta))/Re  # Shear contribution to lift
        '''
        shear_drag_components = ( fields[3, :, t] * 2                    * np.cos(theta) + \
                                 (fields[5, :, t] + fields[4, :, t]) * np.sin(theta))/Re  # Shear contribution to drag
        shear_lift_components = ((fields[5, :, t] + fields[4, :, t]) * np.cos(theta) + \
                                  fields[6, :, t] * 2                    * np.sin(theta))/Re  # Shear contribution to lift


        # Compute the drag and lift components from pressure
        drag_components = -fields[0,:,t] * np.cos(theta)
        lift_components = -fields[0,:,t] * np.sin(theta)

        # Integrate around the circle (sum and multiply by 2 * pi / n for full circle)
        drag[t] = np.sum(drag_components) * dS
        lift[t] = np.sum(lift_components) * dS
        shear_drag[t] = np.sum(shear_drag_components) * dS
        shear_lift[t] = np.sum(shear_lift_components) * dS

    # Total drag and lift (pressure + shear)
    total_drag = drag + shear_drag
    total_lift = lift + shear_lift
    
    results = dict(
        t=times,
        fx_tot=total_drag,
        fx_p=drag,
        fx_visc=shear_drag,
        fy_tot=total_lift,
        fy_p=lift,
        fy_visc=shear_lift,
    )
    # ------------------------------------------------------------------------ #
    # Caching the results
    if cache_dir is not None:
        os.makedirs(os.path.dirname(cache_file), exist_ok=True)
        with open(cache_file, "wb") as f:
            pickle.dump(results, f)

    return results
    
    
    
#
def read_force_torque(file_name: str, cache_dir: str = None) -> dict:
    if cache_dir is not None:
        file_ext = os.path.splitext(file_name)[1]
        cache_file = os.path.join(
            cache_dir, "lift_and_drag",
            os.path.relpath(file_name,
                            cache_dir).replace("../", "").replace("/", "_"))
        cache_file = cache_file.replace(file_ext, ".pkl")

        if os.path.exists(cache_file):
            with open(cache_file, "rb") as f:
                return pickle.load(f)

    # ------------------------------------------------------------------------ #
    # Main computation
    # note!
    # I bet this isn't the best for long logfiles!
    t_list = []
    fx_tot_list = []
    fx_visc_list = []
    fx_p_list = []
    fy_tot_list = []
    fy_visc_list = []
    fy_p_list = []
    
    with open(file_name) as file:
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

    results = dict(
        t = np.array(t_list),
        fx_tot = np.array(fx_tot_list),
        fx_p = np.array(fx_p_list),
        fx_visc = np.array(fx_visc_list),
        fy_tot = np.array(fy_tot_list),
        fy_p = np.array(fy_p_list),
        fy_visc = np.array(fy_visc_list),
        )
    # ------------------------------------------------------------------------ #
    # Caching the results
    if cache_dir is not None:
        os.makedirs(os.path.dirname(cache_file), exist_ok=True)
        with open(cache_file, "wb") as f:
            pickle.dump(results, f)

    return results

def read_brinkman_force(file_name: str, cache_dir: str = None) -> dict:
    if cache_dir is not None:
        file_ext = os.path.splitext(file_name)[1]
        cache_file = os.path.join(
            cache_dir, "lift_and_drag",
            os.path.relpath(file_name,
                            cache_dir).replace("../", "").replace("/", "_"))
        cache_file = cache_file.replace(file_ext, ".pkl")

        if os.path.exists(cache_file):
            with open(cache_file, "rb") as f:
                return pickle.load(f)

    # ------------------------------------------------------------------------ #
    # Main computation
    t_list = []
    fx_list = []
    fy_list = []

    with open(file_name) as file:
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


    results = dict(
        t = np.array(t_list),
        fx_tot = np.array(fx_list),
        fy_tot = np.array(fy_list),
        #note: the pressure viscous split doesn't apply here, we should always
        #check that it exists before we try to plot
        )
    # ------------------------------------------------------------------------ #
    # Caching the results
    if cache_dir is not None:
        os.makedirs(os.path.dirname(cache_file), exist_ok=True)
        with open(cache_file, "wb") as f:
            pickle.dump(results, f)

    return results
