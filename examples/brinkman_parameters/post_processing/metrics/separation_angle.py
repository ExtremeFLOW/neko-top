""""Compute the separation angle of the flow"""

import os
import pickle
import numpy as np
from pynektools.io.read_probes import ProbesReader

# a quick plotting function
def plot_separation_angle(benchmark_results, fig, ax, case, color, linestyle):
    # we only want the largest one
    ax[2].plot(benchmark_results["times"], benchmark_results["angles"][:,0], label =  case["name"], color=color, linestyle=linestyle)
# ============================================================================ #
# Main function for computing the separation angle


def inflection_benchmark(file_name: str, cache_dir: str = None) -> dict:
    """
    Run the inflection benchmark for a given file.

    The function will run the inflection benchmark for the given file and
    return the results.

    Parameters
    ----------
    file_name : str
        The name of the file to run the benchmark on.
    cache_dir : str, optional
        The directory to store the cache files. Cached files will be used if
        they exist. If None, no caching will be used. Default is None.

    Returns
    -------
    benchmark_results : dict
        The results of the benchmark.
        times: np.array
            The time steps.
        angles: np.ndarray
            The angle time series of the inflection points.
        max_freq: np.array
            The dominant frequencies of the inflection points.
        amplitude: np.array
            The amplitudes of the inflection points.
        bias: np.array
            The bias of the inflection points.
        i_boundary: int
            The index at which the boundary regions have formed.
        i_building: int
            The index at which the vortex building is complete.
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
            cache_dir, "inflection_benchmark",
            os.path.relpath(file_name,
                            cache_dir).replace("../", "").replace("/", "_"))

        cache_file = cache_file.replace(file_ext, ".pkl")
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

    # Read in the file and setup the data
    probes = ProbesReader(file_name)

    points = probes.points
    times = probes.times

    fields = np.asarray((probes.fields["u"], probes.fields["v"]))
    del probes

    center = np.array([0.0, 0.0, 0.0])
    angles = track_inflection_point(center, points, fields, times)

    i_boundary, i_building = detect_stable_regions(angles, threshold=1e-6)

    test_times = times[i_building:]
    test_angles = angles[i_building:, :]

    frequencies, spectrum = compute_vortex_frequency(test_times, test_angles)
    max_freq = compute_max_frequency(spectrum, frequencies)
    amplitude = compute_amplitudes(test_angles)
    bias = compute_bias(test_angles)

    benchmark_results = dict(
        times=times,
        angles=angles,
        i_boundary=i_boundary,
        i_building=i_building,
        max_freq=max_freq,
        amplitude=amplitude,
        bias=bias,
    )

    # ------------------------------------------------------------------------ #
    # Caching the results

    if cache_dir is not None:
        os.makedirs(os.path.dirname(cache_file), exist_ok=True)
        with open(cache_file, "wb") as f:
            pickle.dump(benchmark_results, f)

    return benchmark_results


# ============================================================================ #
# Helper functions for computing the separation angle


def compute_bias(angles):
    """
    Compute the bias of the inflection points

    The function will compute the bias of the inflection points. Aka the mean
    value.

    Parameters
    ----------
    angles : np.ndarray
        The angle time series of the inflection points.

    Returns
    -------
    bias : np.array
        The bias of the inflection points.
    """

    bias = np.mean(angles, axis=0)

    return bias


def compute_amplitudes(angles):
    """
    Compute the amplitudes of the inflection points

    The function will compute the amplitudes of the inflection points.

    Parameters
    ----------
    angles : np.ndarray
        The angle time series of the inflection points.

    Returns
    -------
    amplitude : np.array
        The amplitudes of the inflection points.
    """

    mu = np.mean(angles, axis=0)

    amplitude = np.array([
        (angles[:, i] - mu[i]).max() - (angles[:, i] - mu[i]).min()
        for i in range(angles.shape[1])
    ])

    return amplitude


def compute_max_frequency(spectrum, frequencies):
    """
    Compute the dominant frequencies

    The function will compute the dominant frequencies of the inflection points.

    Parameters
    ----------
    spectrum : np.array
        The frequency spectrum of the inflection points.
    frequencies : np.array
        The frequencies list.

    Returns
    -------
    max_freq : np.array
        The dominant frequencies of the inflection points.
    """
    # Print the dominant frequencies
    max_freq = np.zeros(3)
    for i in range(3):
        idx = np.argmax(np.abs(spectrum[:, i]))
        max_freq[i] = np.abs(frequencies[idx])

    return max_freq


def compute_vortex_frequency(times, angles):
    """
    Compute the vortex shedding frequency

    The function will compute the frequency spectrum of the inflection points
    and return the dominant frequencies.

    Parameters
    ----------
    times : np.array
        The time steps.
    angles : np.ndarray
        The angle time series of the inflection points.

    Returns
    -------
    frequencies : np.array
        The frequencies list.
    spectrum : np.array
        The frequency spectrum of the inflection points.
    """

    # Compute the frequency spectrum of the signals
    frequencies = np.fft.fftfreq(len(times), d=times[1] - times[0])
    frequencies = np.fft.fftshift(frequencies)

    spectrum = np.zeros([len(frequencies), 3], dtype=complex)
    for i in range(3):
        spectrum[:, i] = np.fft.fft(angles[:, i] - angles[:, i].mean())
        spectrum[:, i] = np.fft.fftshift(spectrum[:, i])

        # Normalize the spectrum
        spectrum[:, i] /= sum(np.abs(spectrum[:, i]))

    return frequencies, spectrum


def detect_stable_regions(angles, threshold=0.1):
    """
    Detect stable regions in the angle data

    The function will return the indices delimiting the stable regions.

    Parameters
    ----------
    angles : np.ndarray
        The angles of the inflection points.
    threshold : float
        The threshold for the difference between the angles.

    Returns
    -------
    i_boundary : int
        The index at which the boundary regions have formed.
    i_building : int
        The index at which the vortex building is complete.
    """

    # Compute the difference between the inflection points
    diff = np.diff(angles, axis=0)

    i_boundary = 0
    max_diff = 0.0

    for i in range(1, angles.shape[0]):
        if diff[i - 1, :].max() < threshold and max_diff > 1e-2:
            i_boundary = i
            break

        max_diff = max(max_diff, diff[i - 1, :].max())

    i_building = i_boundary
    min_signal = np.array(angles[i_boundary, :], copy=True)
    max_signal = np.array(angles[i_boundary, :], copy=True)
    for i in range(i_boundary, angles.shape[0]):

        criterion_min = min(angles[i, :] - min_signal) < -threshold
        criterion_max = min(max_signal - angles[i, :]) < -threshold

        max_signal = np.maximum(max_signal, angles[i, :])
        min_signal = np.minimum(min_signal, angles[i, :])

        if criterion_min or criterion_max:
            i_building = i

    if i_building == angles.shape[0] - 1:
        i_building = i_boundary

    return i_boundary, i_building


def track_inflection_point(center, points, fields, times):
    """
    Track the inflection point of the flow
    """

    N_points = points.shape[0]
    N_times = times.shape[0]

    angle_inflection = np.zeros([len(times), 3])
    angles = np.zeros(points.shape[0])
    for i_boundary in range(points.shape[0]):
        point = points[i_boundary, :]
        angles[i_boundary] = compute_angle(center, point)

    idx = np.argsort(angles)

    for ti in range(N_times):
        t = times[ti]

        flow = np.zeros(N_points)
        for i_boundary in range(N_points):
            angle = angles[i_boundary]
            field = fields[:, i_boundary, ti]
            flow[i_boundary] = compute_radial_flow(angle, field)

        i0 = int(3 / 4 * len(idx))
        a0 = compute_inflection_angle(i0, 1, flow[idx], angles[idx])

        i0 = int(2 / 4 * len(idx))
        a1 = compute_inflection_angle(i0, -1, flow[idx], angles[idx])

        i0 = int(1 / 4 * len(idx))
        a2 = compute_inflection_angle(i0, 1, flow[idx], angles[idx])

        # The central inflection point cannot be outside the other two.
        if a0 < a1 or a1 < a2:
            a1 = 0.5 * (a0 + a2)

        angle_inflection[ti, :] = [a0, a1, a2]

    return angle_inflection


def compute_inflection_angle(i0, direction, flow, angles) -> float:

    i_p = i0
    i_n = i_p + direction

    while flow[i_p] < 0 or flow[i_n] > 0:
        if flow[i_n] < 0:
            s = -direction
        else:
            s = +direction

        i_p = (i_p + s) % len(angles)
        i_n = (i_n + s) % len(angles)

    # Compute the angle of the inflection point
    angle_p = angles[i_p]
    angle_n = angles[i_n]
    flow_p = flow[i_p]
    flow_n = flow[i_n]

    if abs(flow_n - flow_p) < 1e-8:
        return angle_n
    else:
        return angle_n + (angle_p - angle_n) * flow_n / (flow_n - flow_p)


def compute_radial_flow(angle, field) -> float:
    """
    Compute the radial flow of the field
    """

    # Compute the flow
    flow = -field[0] * np.sin(angle) + field[1] * np.cos(angle)

    return flow


def compute_angle(center, point, axis=np.array([1.0, 0.0, 0.0])):
    """
    Compute the angle between the center and the point
    """

    # Compute the angle
    angle = np.arctan2(point[1] - center[1], point[0] - center[0])

    return angle
