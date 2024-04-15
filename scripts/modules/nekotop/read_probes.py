import csv
import numpy as np


def read_probes(file_name: str) -> tuple:
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

    with open(file_name) as file:
        N_lines = sum(1 for _ in file)

    with open(file_name, "r") as file:
        reader = csv.reader(file)

        # Read in the header
        header = next(reader)
        N_p = int(header[0])
        N_s = int((N_lines - 1 - N_p) / (N_p)) * N_p
        field_names = header[2:]

        points = np.asarray([next(reader) for _ in range(N_p)], dtype=float)
        fields = np.asarray([next(reader) for _ in range(N_s)], dtype=float)

        times = fields[0::N_p, 0]
        fields = fields[:, 1:]

    return (points, fields, times, field_names)
