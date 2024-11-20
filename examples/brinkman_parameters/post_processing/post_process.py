from post_processing_tools import compute_everything, plot_everything, generate_tables
from metrics import *
from functions import *
from experiments import *
import os
import numpy as np

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
    # rings to consider (in the statistics!)
    "force_ring_plot": [0.5, 0.532],
    # rings to consider (in the time series!)
    "force_ring_plot_time": ["051"],
    # rings to plot
    # For taking averages of time series
    "t_start": 100.0  # time to start averaging,
}
# colours
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']


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

i_color = 0
# Loop over each .csv file and run the script for each
for file in experiment_files:
    experiment = experiment_reader(os.path.join(experiments_dir, file))
    experiment_name = file.split('.')[0]

    # ok this will NOT be heavy so we can store them all as a list.
    experiment_tabulated = []

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

        # Construct the path to the case folder
        case_file_path = os.path.join(cases_dir, case["name"])

        # Check if the particular case folder exists
        if os.path.exists(case_file_path):
            print("\tWorking on case:       ", case["name"])
        else:
            print("\tSkipping missing case: ", case["name"])
            continue

        method = case["method"]
        mesh = case["mesh"]
        Re = case["Re"]
        chi = case["chi"]
        radius = case["radius"]

        # use something similar to case for a data point, so we already the
        # names etc
        case_tabulated = case
        data_points = []

        # Force using surface integrals (meshed case)
        color = colors[i_color%len(colors)]
        # -------------------------------------------------------------------- #
        if method == "meshed":
            file_name = os.path.join(case_file_path, 'cylinder.log')
            try:
                force_measure = read_force_torque(file_name, pickle_dir)
                plot_force_measure(force_measure, fig_LD, ax_LD, case, color, '--')
                means, amps = compress_to_data_point(force_measure,
                                                     plot_params)
                data_point = {
                    "mean": means,
                    "amp": amps,
                    "name": "meshed_force"
                }
                data_points.append(data_point)
                del force_measure
            except Exception as e:
                print("Error in case:", case["name"])
                print("\t", e)

        # Force using rho * u (brinkman case)
        # fuck this one...
        # -------------------------------------------------------------------- #
        # if method == "brinkman":
        #     file_name = os.path.join(case_file_path, 'cylinder.log')
        #     try:
        #         force_measure = read_brinkman_force(file_name, pickle_dir)
        #         plot_force_measure(force_measure, fig_LD, ax_LD, case, color, ':')
        #         means, amps = compress_to_data_point(force_measure,
        #                                              plot_params)
        #         data_point = {
        #             "mean": means,
        #             "amp": amps,
        #             "name": "brinkman_force"
        #         }
        #         data_points.append(data_point)
        #         del force_measure
        #     except Exception as e:
        #         print("Error in case:", case["name"])
        #         print("\t", e)

        # Using the rings
        # -------------------------------------------------------------------- #
        for circle in plot_params["force_ring_radii"]:

            try:
                # Lift and Drag
                # ------------------------------------------------------------ #
                file_name = os.path.join(case_file_path,
                                         'circ_' + circle + '.csv')
                force_measure = surface_integral_lift_and_drag(
                    file_name, Re, pickle_dir)

                # append a single curve to the plot
                # just so we don't get a million curves, let's just do the 
                # 050 one.
                if circle in plot_params["force_ring_plot_time"]:
                    plot_force_measure(force_measure, fig_LD, ax_LD, case, color, '-')
                means, amps = compress_to_data_point(force_measure,
                                                     plot_params)
                data_point = {
                    "mean": means,
                    "amp": amps,
                    "name": "force_" + circle
                }
                data_points.append(data_point)
                del force_measure

                # separation angle
                # --------------------------------------------------------#
                sep_angle = inflection_benchmark(file_name, pickle_dir)
                # here we would have some plotting, but for now lets just
                # generate the pickle files
                data_point = {
                    "mean_Str": np.mean(sep_angle['max_freq']),
                    "theta": np.rad2deg(np.max(sep_angle['bias'])),
                    "name": "separation_angle_" + circle
                }
                data_points.append(data_point)

                if circle in plot_params["force_ring_plot_time"]:
                    plot_separation_angle(sep_angle, fig_LD, ax_LD, case, color, '-')
                del sep_angle

            except Exception as e:
                print("Error in case:", case["name"])
                print("\t", e)
            case_tabulated["observables"] = data_points
            experiment_tabulated.append(case_tabulated)
            i_color = i_color + 1

    # here we would finalize all the plots
    os.path.exists(plots_dir) or os.makedirs(plots_dir)
    output_filename = os.path.join(plots_dir,
                                   experiment_name + '_lift_and_drag.png')
    finalize_plot_force_measure(plot_params, fig_LD, ax_LD, output_filename)
    # and we could plot curves based on tabulated data
    if experiment_name == "Filter_radius":
        output_filename = os.path.join(plots_dir, experiment_name)
        plot_study(experiment_tabulated,
                   x_axis_variable="radius",
                   line_axis_variable="method",
                   out_filename=output_filename + '.png')
    if experiment_name == "Re_study":
        output_filename = os.path.join(plots_dir, experiment_name)
        plot_study(experiment_tabulated,
                   x_axis_variable="Re",
                   line_axis_variable="method",
                   out_filename=output_filename + '.png')
    if experiment_name == "Mesh_study":
        output_filename = os.path.join(plots_dir, experiment_name)
        plot_study(experiment_tabulated,
                   x_axis_variable="mesh",
                   line_axis_variable="method",
                   out_filename=output_filename + '.png')
    if experiment_name == "Report_mesh_study_Re200":
        output_filename = os.path.join(plots_dir, experiment_name)
        plot_study(experiment_tabulated,
                   x_axis_variable="mesh",
                   line_axis_variable="method",
                   out_filename=output_filename + '.png')
    if experiment_name == "Report_mesh_study_Re1000":
        output_filename = os.path.join(plots_dir, experiment_name)
        plot_study(experiment_tabulated,
                   x_axis_variable="mesh",
                   line_axis_variable="method",
                   out_filename=output_filename + '.png')

print("All experiments processed.")
