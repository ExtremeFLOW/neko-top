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
        times: np.array
            The time steps.
        points: np.ndarray
            Coordinates of all points on the circle.
        total_drag: np.ndarray
            Total_drag vs time
        drag: np.array
            Pressure drag vs time.
        shear_drag: np.array
            Shear drag vs time.
        total_lift: np.ndarray
            Total_lift vs time
        lift: np.array
            Pressure lift vs time.
        shear_lift: np.array
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

    # Total drag and lift (pressure + shear)(1/Re)(Re=200)
    total_drag = drag + shear_drag
    total_lift = lift + shear_lift
    
    results = dict(
        time_values=times,
        points=points,
        total_drag=total_drag,
        drag=drag,
        shear_drag=shear_drag,
        total_lift=total_lift,
        lift=lift,
        shear_lift=shear_lift,
    )
    # ------------------------------------------------------------------------ #
    # Caching the results
    if cache_dir is not None:
        os.makedirs(os.path.dirname(cache_file), exist_ok=True)
        with open(cache_file, "wb") as f:
            pickle.dump(results, f)

    return results
    
    
    
#
