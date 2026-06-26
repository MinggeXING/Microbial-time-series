#!/bin/bash
# run_edm_am_cd0.sh - Run MDR CV tasks in Conda r433 environment

# Load Conda initialization
source /home/xmg/miniconda3/enter/etc/profile.d/conda.sh

# Activate the specific environment
conda activate r433

# Check if environment activation succeeded
if [[ -z "$CONDA_PREFIX" ]]; then
    echo "$(date) - Error: Failed to activate Conda environment 'r433'"
    exit 1
fi

# Get the R interpreter path from the environment
R_PATH="$CONDA_PREFIX/bin/R"

# Check if R exists
if [[ ! -f "$R_PATH" ]]; then
    echo "$(date) - Error: R interpreter not found in the environment"
    exit 1
fi

# Change to working directory
cd /home/xmg/work || exit 1

# Run R script and log output
echo "$(date) - Starting MDR CV task"
nohup "$R_PATH" --vanilla < EDM_AM_Cd0.R >> edm_am_cd0.log 2>&1 &