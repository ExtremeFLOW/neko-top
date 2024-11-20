import os

# Local imports
from metrics import *
from functions.readers import experiment_reader
from experiments import *

# Define useful paths for the script
current_dir = os.path.dirname(os.path.abspath(__file__))
root_dir = os.path.abspath(os.path.join(current_dir, "../../.."))
experiments_dir = os.path.join(root_dir, "examples", "brinkman_parameters",
                               "experiments")
cases_dir = os.path.join(root_dir, "results", "brinkman_parameters", "cases")
#cases_dir = os.path.join(root_dir, "results", "only_stats","brinkman_parameters", "cases")

# Create logs and plots folders in the experiments directory if they don't exist
plots_dir = os.path.abspath(os.path.join(current_dir, "../", "plots"))
tables_dir = os.path.abspath(os.path.join(current_dir, "../", "tables"))
cache_dir = os.path.join(root_dir, "results", "brinkman_parameters", "cache")
#cache_dir = os.path.join(root_dir, "results", "only_stats", "brinkman_parameters", "cache")

# ============================================================================ #
# Post process the experiments

plot_params = {
    # Settings for post processing
    # Lift and drag
    "if_lift_and_drag": True,
    "lift_axis": [-1, 1],  # axis for plotted lift,
    "drag_axis": [0, 1.5],  # axis for plotted drag,
    # Statistics interpolation (wake deficit)
    "if_stats_wake": True,
    "wake_y_lim": [-15, 15],  # y limits for wake profiles,
    "wake_n_pts": 300,  # number points in the wake,
    "wake_positions": [10, 15, 20],  # position for wakes,
    "wake_U_lim": [0.5, 1.2],  # axis limits for wakes plot,
    # Statistics interpolation (force rings)
    "if_stats_force": True,
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

# List all .csv files in the experiments directory
experiment_files = [
    f for f in os.listdir(experiments_dir) if f.endswith('.csv')
]

# Loop over each .csv file and run the script for each
for file in experiment_files:
    i_color = 0
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
    # example: plot_lift_and_drag(experiment, 'Re', cache_dir)
    #
    # Possibly which should also accept a list of plot options so we can control
    # styling etc.

    # set up the axis for a test plot
    fig_LD, ax_LD = init_plot_force_measure(plot_params)
    fig_wake, ax_wake = init_plot_a_wake_line(plot_params)
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

        # somehow I've still fucked up this color order...
        if i_color == 0:
            ind = 0
        else:
            ind = i_color % len(colors)
        color = colors[ind]
        # use something similar to case for a data point, so we already the
        # names etc
        case_tabulated = case
        data_points = []

        # Wake lines
        # -------------------------------------------------------------------- #
        try:
            wake_line = wake(case_file_path, plot_params, cache_dir)
            wake_line["name"] = case["name"]
            plot_a_wake_line(wake_line, plot_params, ax_wake, fig_wake, color, '-')
        except Exception as e:
            print("Error in statistics case:", case["name"])
            print("\t", e)
        # -------------------------------------------------------------------- #
        # for processing the stats files, I guess we could also use the cache.


        # Force using surface integrals (meshed case)
        # -------------------------------------------------------------------- #
        if method == "meshed":
            file_name = os.path.join(case_file_path, 'cylinder.log')
            try:
                force_measure = read_force_torque(file_name, cache_dir)
                plot_force_measure(force_measure, fig_LD, ax_LD, case, color,
                                   '--')
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
        #         force_measure = read_brinkman_force(file_name, cache_dir)
        #         plot_force_measure(force_measure, fig_LD, ax_LD)
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
                    file_name, Re, cache_dir)

                # append a single curve to the plot
                # just so we don't get a million curves, let's just do the
                # 050 one.
                if circle in plot_params["force_ring_plot_time"]:
                    plot_force_measure(force_measure, fig_LD, ax_LD, case,
                                       color, '-')
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
                sep_angle = inflection_benchmark(file_name, cache_dir)
                # here we would have some plotting, but for now lets just
                # generate the pickle files
                data_point = {
                    "mean_Str": np.mean(sep_angle['max_freq']),
                    "theta": np.rad2deg(np.max(sep_angle['bias'])),
                    "name": "separation_angle_" + circle
                }
                data_points.append(data_point)

                if circle in plot_params["force_ring_plot_time"]:
                    plot_separation_angle(sep_angle, fig_LD, ax_LD, case,
                                          color, '-')
                del sep_angle

            except Exception as e:
                print("Error in case:", case["name"])
                print("\t", e)
        # -------------------------------------------------------------------- #
        # tabulate them all and move on
        case_tabulated["observables"] = data_points
        experiment_tabulated.append(case_tabulated)
        i_color = i_color + 1
        # -------------------------------------------------------------------- #
        # for processing the stats files, I guess we could also use the cache.

    # here we would finalize all the plots
    os.path.exists(plots_dir) or os.makedirs(plots_dir)
    output_filename = os.path.join(plots_dir,
                                   experiment_name + '_lift_and_drag.png')
    finalize_plot_force_measure(plot_params, fig_LD, ax_LD, output_filename)
    os.path.exists(plots_dir) or os.makedirs(plots_dir)
    output_filename = os.path.join(plots_dir,
                                   experiment_name + '_wake_line.png')
    finalize_plot_a_wake_line(plot_params, ax_wake, fig_wake, output_filename)
    # and we could plot curves based on tabulated data
    if experiment_name == "Filter_radius":
        output_filename = os.path.join(plots_dir, experiment_name)
        ax, fig = plot_study(experiment_tabulated,
                   x_axis_variable="radius",
                   line_axis_variable="method")
        # extract just the meshed case
        case_file_path = os.path.join(cases_dir, "meshed_mesh_4_re_200")
        # extract the forces
        file_name = os.path.join(case_file_path, 'cylinder.log')
        force_measure = read_force_torque(file_name, cache_dir)
        means, amps = compress_to_data_point(force_measure, plot_params)
        file_name = os.path.join(case_file_path, 'circ_' + '050' + '.csv')
        sep_angle = inflection_benchmark(file_name, cache_dir)
        mean_Str = np.mean(sep_angle['max_freq'])
        theta = np.rad2deg(np.max(sep_angle['bias']))
        ax[0].plot([0.01,0.1],[means["fx_tot"],means["fx_tot"]], 
            linestyle = '--', color = 'k', label = "meshed reference")
        ax[0].legend() 
        ax[1].plot([0.01,0.1],[theta, theta], linestyle = '--', color = 'k', 
            label = "meshed reference")
        ax[1].legend() 
        ax[2].plot([0.01,0.1],[mean_Str, mean_Str], linestyle = '--', 
            color = 'k', label = "meshed reference")
        ax[2].legend() 
        fig.savefig(output_filename)
    if experiment_name == "Re_study":
        output_filename = os.path.join(plots_dir, experiment_name)
        ax, fig =  plot_study(experiment_tabulated,
                   x_axis_variable="Re",
                   line_axis_variable="method")
        fig.savefig(output_filename)
    if experiment_name == "Mesh_study":
        output_filename = os.path.join(plots_dir, experiment_name)
        ax, fig = plot_study(experiment_tabulated,
                   x_axis_variable="mesh",
                   line_axis_variable="method")
        fig.savefig(output_filename)
    if experiment_name == "Report_mesh_study_Re200":
        output_filename = os.path.join(plots_dir, experiment_name)
        ax, fig = plot_study(experiment_tabulated,
                   x_axis_variable="mesh",
                   line_axis_variable="method")
        fig.savefig(output_filename)
    if experiment_name == "Report_mesh_study_Re1000":
        output_filename = os.path.join(plots_dir, experiment_name)
        ax, fig = plot_study(experiment_tabulated,
                   x_axis_variable="mesh",
                   line_axis_variable="method")
        fig.savefig(output_filename)

print("All experiments processed.")
