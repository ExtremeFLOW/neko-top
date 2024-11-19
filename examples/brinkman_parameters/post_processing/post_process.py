from post_processing_tools import compute_everything, plot_everything, generate_tables
from metrics import *
from functions import *
import os
import csv
import matplotlib.pyplot as plt

# OK, I (Harry) am definitely not the person to be making a judgement call like
# this hahahaha, it's very computer science-y and is based on a conversation
# over coffee with Bertie and a very superficial amount of googling.
#
# The current csv format is always going to be a killer for us to read, and if
# we have to read them all, every time, just to fuck around with the axis of
# a plot it's going to be a nightmare.
#
# I'm SURE there's a better way of doing this, but according to google, Pickle
# is a good option.
#
# So we divide this script into a preprocessing reader + pickle saver
#
# And then we have the actual plotting script.
import pickle

plot_params = {
    # Settings for post processing
    # Lift and drag
    "if_lift_and_drag": True,
    "lift_axis": [-1, 1],  # axis for plotted lift,
    "drag_axis": [0, 1.5],  # axis for plotted drag,
    # Statistics interpolation (wake deficit)
    "if_stats_wake": False,
    "wake_y_lim": [-15, 15],  # y limits for wake profiles,
    "wake_n_pts": 300,  # number points in the wake,
    "wake_positions": [10, 15, 20],  # position for wakes,
    "wake_U_lims": [0.5, 1.2],  # axis limits for wakes plot,
    # Statistics interpolation (force rings)
    "if_stats_force": False,
    "force_n_pts": 360,  # number of interpolation points,
    "force_ring_radii": ["050", "0501", "0502", "0505", "051", "052", "055"],
    # rings to consider
    "force_ring_plot": [0.5, 0.532],
    # rings to plot
    # For taking averages of time series
    "t_start": 100.0  # time to start averaging,
}

# Define useful paths for the script
current_dir = os.path.dirname(os.path.abspath(__file__))
root_dir = os.path.abspath(os.path.join(current_dir, "../../.."))
experiments_dir = os.path.join(root_dir, "examples", "brinkman_parameters",
                               "experiments")
cases_dir = os.path.join(root_dir, "results", "brinkman_parameters", "cases")

# Create logs and plots folders in the experiments directory if they don't exist
plots_dir = os.path.abspath(os.path.join(current_dir, "../", "plots"))
tables_dir = os.path.abspath(os.path.join(current_dir, "../", "tables"))
pickle_dir = os.path.join(root_dir, "results", "brinkman_parameters", "cache")

# List all .csv files in the experiments directory
experiment_files = [
    f for f in os.listdir(experiments_dir) if f.endswith('.csv')
]

# Loop over each .csv file and run the script for each
for file in experiment_files:
    experiment = experiment_reader(os.path.join(experiments_dir, file))
    experiment_name = file.split('.')[0]

    print("Starting with the experiment:", experiment_name)

    # This should be encased in a function.
    #
    # This function should take the experiment list, the variable we want to
    # plot against and the pickle directory as arguments.
    #
    # example: plot_lift_and_drag(experiment, 'Re', pickle_dir)
    #
    # Possibly which should also accept a list of plot options so we can control
    # styling etc.

    # set up the axis for a test plot
    fig_LD, ax_LD = init_plot_force_measure(plot_params)
    # Loop through each case/case in the experiment
    for case_number, case in enumerate(experiment):
        print("Working on case number:", case_number, "   ", case["name"])

        # Construct the path to the case folder
        case_file_path = os.path.join(cases_dir, case["name"])

        # Check if the particular case folder exists
        if os.path.exists(case_file_path):
            method = case["method"]
            mesh = case["mesh"]
            Re = case["Re"]
            chi = case["chi"]
            radius = case["radius"]

            for circle in plot_params["force_ring_radii"]:

                try:
                    # either Calculate or load the case and append case data
                    # this will be full of all the functions we're interested in,
                    # for now it's just one to test
                    file_name = os.path.join(case_file_path,
                                             'circ_' + circle + '.csv')
                    force_measure = surface_integral_lift_and_drag(
                        file_name, Re, pickle_dir)

                    # append a single curve to the plot
                    plot_force_measure(force_measure, fig_LD, ax_LD)
                    del force_measure
                except Exception as e:
                    print("Error in case:", case["name"])
                    print("\t", e)

    # here we would finalize all the plots
    output_filename = os.path.join(plots_dir,
                                   experiment_name + '_lift_and_drag.png')
    finalize_plot_force_measure(plot_params, fig_LD, ax_LD, output_filename)

print("All experiments processed.")
