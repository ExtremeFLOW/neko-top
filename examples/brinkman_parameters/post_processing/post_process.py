from post_processing_tools import compute_everything, plot_everything, generate_tables
from metrics import *
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
plot_params = {}
plot_params["precompute"]       = True       # do we read the csv or the pickle?

# Settings for post processing
# Lift and drag
plot_params["if_lift_and_drag"] = True
plot_params["lift_axis"]        = [-1, 1]    # axis for plotted lift   
plot_params["drag_axis"]        = [0, 1.5]   # axis for plotted drag   
# Statistics interpolation (wake deficit)
plot_params["if_stats_wake"]    = False
plot_params["wake_y_lim"]       = [-15, 15]  # y limits for wake profiles
plot_params["wake_n_pts"]       = 300        # number points in the wake
plot_params["wake_positions"]   = [10,15,20] # position for wakes
plot_params["wake_U_lims"]      = [0.5, 1.2] # axis limits for wakes plot
# Statistics interpolation (force rings)
plot_params["if_stats_force"]   = False
plot_params["force_n_pts"]      = 360        # number of interpolation points
plot_params["force_ring_radii"] = [0.50, 0.501, 0.502, 0.504, 0.508, 0.516, 0.532]
                                             # rings to consider
plot_params["force_ring_plot"]  = [0.5, 0.532]
                                             # rings to plot
# For taking averages of time series
plot_params["t_start"]          = 100.0      # time to start averaging


# Get the current directory
current_dir = os.path.dirname(os.path.abspath(__file__))

# Navigate to the parent directory where the experiments and cases folders are located
parent_dir = os.path.abspath(os.path.join(current_dir, "../../.."))

# Define the path to 'experiments'
experiments_dir = os.path.join(parent_dir, "results", "brinkman_parameters_fake", "experiments")
# experiments_dir = os.path.join(parent_dir, "results", "only_stats", "brinkman_parameters", "experiments")

# Define the path to 'cases'
cases_dir = os.path.join(parent_dir, "results", "brinkman_parameters_fake", "cases")
# cases_dir = os.path.join(parent_dir, "results", "only_stats", "brinkman_parameters", "cases")

# Create logs and plots folders in the experiments directory if they don't exist
logs_dir = os.path.join(experiments_dir, "logs")
plots_dir = os.path.join(experiments_dir, "plots")
tables_dir = os.path.join(experiments_dir, "tables")
pickle_dir = os.path.join(experiments_dir, "pickle_files")
os.makedirs(logs_dir, exist_ok=True)
os.makedirs(plots_dir, exist_ok=True)
os.makedirs(tables_dir, exist_ok=True)
os.makedirs(pickle_dir, exist_ok=True)



# List all .csv files in the experiments directory
csv_files = [f for f in os.listdir(experiments_dir) if f.endswith('.csv')]

# Loop over each .csv file and run the script for each
for csv_file in csv_files:
    experiment = os.path.splitext(csv_file)[0]  # Remove the .csv extension
    file_path = os.path.join(experiments_dir, csv_file)
    output_file_path = os.path.join(logs_dir, experiment + "_log.csv")

    # set up the axis for a test plot
    fig_LD, ax_LD = init_plot_force_measure(plot_params)

    print("Starting with the experiment:", experiment)
    if os.path.exists(file_path):
        # Open the original CSV file and read it
        with open(file_path, 'r') as file:
            reader = csv.DictReader(file)  # Read the file as a dictionary (each row is a dictionary)

            # Get the fieldnames (headers) from the original CSV
            fieldnames = reader.fieldnames + ['status']  # Add 'status' to the fieldnames
            
            # Open the output CSV file for writing
            with open(output_file_path, 'w', newline='') as output_file:
                writer = csv.DictWriter(output_file, fieldnames=fieldnames)
                writer.writeheader()  # Write the header row

                case_number = 1
                # Loop through each row/case in the experiment
                for row in reader:
                    print("Working on case number:", case_number, "   ", row[reader.fieldnames[0]])
                
                    # Get the folder name from the first column (this is the experiment case folder)
                    case_folder_name = row[reader.fieldnames[0]]

                    # Construct the path to the case folder
                    case_file_path = os.path.join(cases_dir, case_folder_name)

                    # Check if the particular case folder exists
                    if os.path.exists(case_file_path):
                        row['status'] = 'exists'  # If folder exists, set status as 'exist'
                        method = row[reader.fieldnames[1]]
                        mesh = int(row[reader.fieldnames[2]])
                        Re = float(row[reader.fieldnames[3]])
                        if method == 'meshed':
                            #If method is meshed then we can't read the following values from the file
                            # therefore, I am just giving them some extream initial values
                            chi = -1.0
                            implicit = False
                            radius = -1.0
                        else:
                            chi = float(row[reader.fieldnames[4]])
                            implicit = bool(row[reader.fieldnames[5]])
                            radius = float(row[reader.fieldnames[6]])
                        
                        # either Calculate or load the case and append case data
                        # this will be full of all the functions we're interested in,
                        # for now it's just one to test
                        file_name = os.path.join(case_file_path, 'circ_0501.csv')
                        cache_file = pickle_dir
                        force_measure = surface_integral_lift_and_drag(file_name, Re, cache_file)
                        # append a single curve to the plot
                        plot_force_measure(force_measure, fig_LD, ax_LD)
                        del force_measure
                        
                    else:
                        row['status'] = 'not exists'  # If folder does not exist, set status as 'not exist'

                    # Write the row with the new 'status' field to the output CSV
                    writer.writerow(row)
                    case_number += 1

    else:
        print(f"The file {file_path} does not exist.")
    # here we would finalize all the plots
    output_filename = plots_dir + '/' + experiment + '_lift_and_drag.png'
    finalize_plot_force_measure(plot_params, fig_LD, ax_LD, output_filename)

print("All experiments processed.")
