#!/bin/bash
#SBATCH --job-name=oyster_qmra
#SBATCH --output=oyster_qmra_%j.log
#SBATCH --error=oyster_qmra_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=weifan.wu@ufl.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128gb
#SBATCH --time=12:00:00
#SBATCH --account=dianataft
#SBATCH --qos=dianataft

# Load R 4.3 (contains user's R packages: tidyverse, mc2d, brms)
module purge
module load R/4.3

# Set environment variables
export R_MAX_VSIZE=100Gb
export FORCE_RECOMPUTE=TRUE

echo "Starting Oyster Norovirus QMRA run on HyperGator..."
echo "Job ID: $SLURM_JOB_ID"
echo "Running on node: $SLURM_NODELIST"
date

# Render/run the Quarto document (executes full QMRA pipeline and saves all RDS/CSV outputs to results/)
Rscript -e "options(bitmapType='cairo'); knitr::purl('scripts/oyster_qmra_master.qmd', output='scratch_run.R'); source('scratch_run.R')"

echo "QMRA job completed successfully!"
date
