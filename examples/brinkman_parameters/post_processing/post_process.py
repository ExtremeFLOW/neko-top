from post_processing_tools import compute_everything, plot_everything
import os
import csv
import matplotlib.pyplot as plt



experiment = "Filter_radius"

# Get the current directory
current_dir = os.path.dirname(os.path.abspath(__file__))

# Navigate to the 'cases' folder by going back and then into 'results/brinkman_parameters/cases'
parent_dir = os.path.abspath(os.path.join(current_dir, "../../.."))

# Define the path to 'experiments'
# experiments_dir = os.path.join(parent_dir, "results", "brinkman_parameters_fake", "experiments")
experiments_dir = os.path.join(parent_dir, "results", "only_stats", "brinkman_parameters", "experiments")

# Define the path to 'cases'
cases_dir = os.path.join(parent_dir, "results", "only_stats", "brinkman_parameters", "cases")

# Path to the original file
file_path = os.path.join(experiments_dir, experiment + ".csv")

# Path to the output file
output_file_path = os.path.join(experiments_dir, experiment + "_log.csv")



#results_path = "../../../results/brinkman_parameters/cases/"
#meshed_path = "meshed_mesh_2_re_2000"
#brinkman_path = "brinkman_mesh_2_re_200_chi_100_radius_0-03"
#meshed = {}
#IBM = {}
#meshed["name"] = 'meshed'
#IBM["name"] = 'IBM'


case_list = []



print("Starting with the experiment ", experiment)
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
            # Loop through each row/cases in the experiment
            for row in reader:
                print("Working on case number:" , case_number, "   ", row[reader.fieldnames[0]] )
            
                # Get the folder name from the first column (this is the experiment case folder)
                case_folder_name = row[reader.fieldnames[0]]

                # Construct the path to the case folder
                case_file_path = os.path.join(cases_dir, case_folder_name)

                # Check if the particular case folder exists
                if os.path.exists(case_file_path):
                    row['status'] = 'exist'  # If folder exists, set status as 'exist'
                    method = row[reader.fieldnames[1]]
                    mesh = int(row[reader.fieldnames[2]])
                    Re = float(row[reader.fieldnames[3]])
                    chi = float(row[reader.fieldnames[4]])
                    implicit = bool(row[reader.fieldnames[5]])
                    radius = float(row[reader.fieldnames[6]])
                    # these ones are strange, they belong to Niclas's stuff
                    # rmax = float(row[reader.fieldnames[7]])
                    # rpower = float(row[reader.fieldnames[8]])
                    
                    case = compute_everything(cases_dir + "/", case_folder_name, method, mesh, Re, chi, implicit, radius)
                    case_list.append(case)
                    
                else:
                    row['status'] = 'not exist'  # If folder does not exist, set status as 'not exist'

                # Write the row with the new 'status' field to the output CSV
                writer.writerow(row)
                case_number = case_number+1

    #print(f"New file with status column has been created: {output_file_path}")

else:
    print(f"The file {file_path} does not exist.")

plot_everything(case_list, experiments_dir, experiment)
