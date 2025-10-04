#!/bin/bash

#SBATCH --partition=sapphire,intermediate,hsph
#SBATCH --qos=normal
#SBATCH --account=dominici_lab
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --job-name baer
#SBATCH --mail-type ALL             
#SBATCH --mail-user tbui@hsph.harvard.edu
#SBATCH --output fasrc_run.out
#SBATCH --mem=16GB
#SBATCH --time=24:00:00
#SBATCH --array=1-1000

# set R packages and rstudio server singularity image locations
my_packages=${HOME}/R/ifxrstudio/RELEASE_3_18
rstudio_singularity_image="/n/singularity_images/informatics/ifxrstudio/ifxrstudio:RELEASE_3_18.sif"

export OMP_NUM_THREADS=1
export SCENARIO=$1

# run myscript.R using RStudio Server signularity image
singularity exec --env R_LIBS_USER=${my_packages} ${rstudio_singularity_image} Rscript '2. sim run.R' > logs/sim_scen$1_array${SLURM_ARRAY_TASK_ID}.Rout 2>&1
