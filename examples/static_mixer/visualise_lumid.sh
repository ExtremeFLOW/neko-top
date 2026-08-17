#!/bin/bash

# In this file make changes to the SBATCH variables to control the LUMI
# hpc settings.
#
# In addition, all modules should be loaded and python virtualenv should be
# setup if python is used in either testing or visualisation.
# Modules and python setups can be done in a separate file and supplied through
# the FILES variable in submit.sh. This will ensure a uniform setup.

# =============================================================================
# Define the SBATCH options here.

# --  Technical Options

# Queue name
#SBATCH --job-name=static_mixer_vis
#SBATCH --account=project_465002485
#SBATCH --partition=lumid
#SBATCH --dependency=singleton

# Ask for n cores placed on R host.
#SBATCH --tasks=8
#SBATCH --cpus-per-task=1
#SBATCH --gpus=1
#SBATCH --mem=48G

# Time specifications (dd-hh:mm:ss)
#SBATCH --time 00-04:00:00

# -- Notification options

# Set the email to receive to and when to receive it
#SBATCH --mail-type=END    # Send notification at completion

# -- Mandatory options, change with great care.

# Definitions of output files.
#SBATCH --output vis_out.log
#SBATCH --error vis_err.log

# ============================================================================ #
# Determine if the script is run on the HPC or locally

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

source ~/Software/load.sh
cd ~/Projects/neko-top

result_root=results/static_mixer
gifmaker=scripts/python/gif_maker.py
visualiser=examples/static_mixer/visualise.py

experiments=(
    petsc/128_re_1000_pe_1000
    petsc/128_re_1500_pe_1500
    petsc/64_re_1000_pe_1000
    petsc/128_re_3000_pe_3000
    petsc/64_re_1500_pe_1500
    petsc/64_re_3000_pe_3000
    petsc_phmg/128_re_1000_pe_1000
    petsc_phmg/128_re_1500_pe_1500
    petsc_phmg/128_re_3000_pe_3000
    petsc_phmg/64_re_1000_pe_1000
    petsc_phmg/64_re_1500_pe_1500
    petsc_phmg/64_re_3000_pe_3000
)
colors=(
    black
    white
)

for exp in "${experiments[@]}"; do
    if [[ ! -d "$result_root/$exp" ]]; then
        echo "Results for experiment $exp not found. Skipping visualization."
        continue
    else
        echo "Visualizing results for experiment $exp."
    fi

    for color in "${colors[@]}"; do
        Temp_gif=$result_root/$exp/visualisation/Temperature_${color}.gif
        Vel_gif=$result_root/$exp/visualisation/Velocity_${color}.gif
        if [[ ! -f "$Temp_gif" || ! -f "$Vel_gif" ]]; then

            # Run the ParaView visualizer.
            srun -u pvbatch $visualiser --text-color=$color $result_root/$exp/

            # If the srun was cancelled, skip the gif creation to avoid errors
            if [[ $? -ne 0 ]]; then
                echo "Visualization for $exp with color $color was cancelled."
                continue
            fi

            # Create the gifs from the generated png files.
            python $gifmaker --pattern Temperature_[0-9]*.png --fps 25 \
                --quality high $result_root/$exp/visualisation/${color}/
            python $gifmaker --pattern Velocity_[0-9]*.png --fps 25 \
                --quality high $result_root/$exp/visualisation/${color}/
        fi
    done
done


# ==============================   End of File   ==============================
