import csv
import os


def experiment_reader(file_name: str):
    """
    Opens ad reads the CSV based experiment files.
    The CSV file must have the following fields:
    name: The name of the experiment
    method: The method used in the experiment
    mesh: The mesh size used in the experiment
    The other fields are the parameters of the experiment and should be float
    values.
    Missing values are represented by '-' and are replaced by None.

    Parameters
    ----------
    file_name : str
        The path to the CSV file.

    Returns
    -------
    list
        A list of dictionaries with the experiment data.
        Each entry in the list is a dictionary with the fields of the CSV file.
    """

    # Open the file
    with open(file_name, 'r') as file:
        reader = csv.reader(file)

        # The header defines the fields of the CSV file
        header = next(reader)

        experiment = []
        # Read the data
        for row in reader:

            # Build the dictionary
            ex = {}
            for i, field in enumerate(header):
                if field == 'name':
                    ex[field] = row[i]
                elif field == 'method':
                    ex[field] = row[i]
                elif field == 'mesh':
                    ex[field] = int(row[i])
                else:
                    # Missing values '-' are replaced by None
                    ex[field] = float(row[i]) if row[i] != '-' else None

            experiment.append(ex)

    return experiment
