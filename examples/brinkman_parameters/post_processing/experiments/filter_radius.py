# Update paths to import the functions

# ============================================================================ #
# Import packages
import matplotlib.pyplot as plt
import os
import sys

try:
    from functions.experiment_reader import experiment_reader
except ImportError or ModuleNotFoundError:
    sys.path.append(os.path.join(os.path.dirname(__file__), '../functions'))
    from experiment_reader import experiment_reader
try:
    from metrics.separation_angle import inflection_benchmark
except ImportError or ModuleNotFoundError:
    sys.path.append(os.path.join(os.path.dirname(__file__), '../metrics'))
    from separation_angle import inflection_benchmark

# ============================================================================ #


def filter_radius(experiment_folder: str, results_folder: str):

    cache_dir = os.path.join(results_folder, 'cache')

    # Read the experiment file
    experiments = experiment_reader(
        os.path.join(experiment_folder, 'Filter_radius.csv'))

    brinkman_data = {}
    idw_data = {}

    for experiment in experiments:
        case_folder = os.path.join(results_folder, experiment['name'])
        circ_051 = os.path.join(case_folder, 'circ_051.csv')

        if not os.path.exists(circ_051):
            print(
                f"File {circ_051} not found. Skipping {experiment['name']}\n")
            continue

        if experiment["method"] == "brinkman":

            separation_angles = inflection_benchmark(circ_051, cache_dir)

            # Append the data to the plot
            brinkman_data[
                experiment['radius']] = separation_angles['max_freq'][0]

        if experiment["method"] == "idw":
            separation_angles = inflection_benchmark(circ_051, cache_dir)

            # Append the data to the plot
            idw_data[experiment['rmax']] = separation_angles['max_freq'][0]

    # Add the legend
    # Create a figure and axis
    fig, ax = plt.subplots()

    # Plot the data
    ax.plot(list(brinkman_data.keys()),
            list(brinkman_data.values()),
            label='Brinkman')
    ax.plot(list(idw_data.keys()), list(idw_data.values()), label='IDW')

    ax.legend()
    plt.show(block=True)


if __name__ == '__main__':

    # Get current directory
    current_dir = os.path.dirname(os.path.realpath(__file__))
    root_folder = os.path.join(current_dir, '../../')
    experiment_folder = os.path.join(root_folder, 'experiments')
    results_folder = os.path.join(
        current_dir, '../../../../results/hpc/brinkman_parameters/cases/')

    experiment_folder = os.path.abspath(experiment_folder)
    results_folder = os.path.abspath(results_folder)

    filter_radius(experiment_folder, results_folder)
