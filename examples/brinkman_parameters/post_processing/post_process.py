from post_processing_tools import compute_everything, plot_everything, generate_tables
import os
import csv
import matplotlib.pyplot as plt


# Settings for post processing
plot_params = {}
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
experiments_dir = os.path.join(parent_dir, "results", "brinkman_parameters", "experiments")
# experiments_dir = os.path.join(parent_dir, "results", "only_stats", "brinkman_parameters", "experiments")

# Define the path to 'cases'
cases_dir = os.path.join(parent_dir, "results", "brinkman_parameters", "cases")
# cases_dir = os.path.join(parent_dir, "results", "only_stats", "brinkman_parameters", "cases")

# Create logs and plots folders in the experiments directory if they don't exist
logs_dir = os.path.join(experiments_dir, "logs")
plots_dir = os.path.join(experiments_dir, "plots")
tables_dir = os.path.join(experiments_dir, "tables")
os.makedirs(logs_dir, exist_ok=True)
os.makedirs(plots_dir, exist_ok=True)
os.makedirs(tables_dir, exist_ok=True)

# List all .csv files in the experiments directory
csv_files = [f for f in os.listdir(experiments_dir) if f.endswith('.csv')]

# Loop over each .csv file and run the script for each
for csv_file in csv_files:
    experiment = os.path.splitext(csv_file)[0]  # Remove the .csv extension
    file_path = os.path.join(experiments_dir, csv_file)
    output_file_path = os.path.join(logs_dir, experiment + "_log.csv")

    case_list = []

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
                        
                        # Calculate and append case data
                        case = compute_everything(cases_dir + "/", case_folder_name, method, mesh, Re, chi, implicit, radius, plot_params)
                        case_list.append(case)
                        
                    else:
                        row['status'] = 'not exists'  # If folder does not exist, set status as 'not exist'

                    # Write the row with the new 'status' field to the output CSV
                    writer.writerow(row)
                    case_number += 1

    else:
        print(f"The file {file_path} does not exist.")
    plot_everything(case_list, plots_dir, experiment, plot_params)
    generate_tables(case_list, tables_dir, experiment, plot_params)
print("All experiments processed.")

