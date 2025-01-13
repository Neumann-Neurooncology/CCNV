#!/bin/bash
#------------------------------------------------------
#SBATCH --job-name=capper_ccnv                      ### the name of your job (and maybe put your name in there as well)
#SBATCH --output=%x_%A_%a.out                           ### the output file for stdout (see: `man sbatch` for format) (will be placed in --chdir or on absolute path)
#------------------------------------------------------
#SBATCH --partition=small                        ### request partition
#SBATCH --ntasks=1                                      ### request number of tasks per node
#SBATCH --cpus-per-task=4 		### number of threads per task (OMP threads)
#------------------------------------------------------
##SBATCH --mail-type=ARRAY_TASKS,END,FAIL                ### some job notifications to receive (NONE, BEGIN, END, FAIL, REQUEUE, ALL), ARRAY_TASKS will send the mail on each task
#------------------------------------------------------
#SBATCH --time=0-72:00:00                               ### total run time limit
#------------------------------------------------------
#SBATCH --array=1-92%30                                  ### specify the integer range and the number of parallel running array jobs (limited to 1000 on HSUper) (it is better to cluster tasks to one array task if they are very small)

# it is a good practice to kill the script if anything weird happens
set -e

echo "#=============================================================================#"
# load python together with numpy via the module environment system
# currently, py-numpy is only available compiled with gcc@12.1.0
# hence, gcc must be loaded beforehand
# this might change with future HSUper software updates
cd "/beegfs/home/s/schumany/AntoniaCCNV"
echo "#=============================================================================#"
echo "    job started: " $(date)
echo "#=============================================================================#"
echo "           user: " $USER
echo "         run at: " $PWD
echo "   array job ID: " $SLURM_ARRAY_JOB_ID
echo "  array TASK ID: " $SLURM_ARRAY_TASK_ID
echo "         job ID: " $SLURM_JOB_ID
echo "       job name: " $SLURM_JOB_NAME
echo "      partition: " $SLURM_JOB_PARTITION
echo "#=============================================================================#"

apptainer exec ubuntu_ccnv.sif Rscript /beegfs/home/s/schumany/AntoniaCCNV/Scripts/capper_template.R $SLURM_ARRAY_TASK_ID

