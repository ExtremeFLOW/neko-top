import csv
import numpy as np


def read_probes(file_name) -> tuple:
    """
    Read in the probes file and return the points, fields, times and field
    names.

    Parameters
    ----------
    file_name : str
        The name of the file to read in.

    Returns
    -------
    points : np.array
        The points in the domain.
    fields : np.array
        The fields at each time step.
    times : np.array
        The times at each time step.
    field_names : list
        The names of the fields.
    """

    with open(file_name) as f:
        N_lines = sum(1 for _ in f)

    with open(file_name, 'r') as f:
        reader = csv.reader(f)

        # Read in the header
        header = next(reader)
        N_points = int(header[0])
        field_names = header[2:]

        N_times = (N_lines - 1 - N_points)

        points = np.array([next(reader) for _ in range(N_points)], dtype=float)
        fields = np.array([next(reader) for _ in range(N_times)], dtype=float)
        times = fields[1::N_points, 0]
        fields = fields[:, 1:]

    return points, fields, times, field_names
