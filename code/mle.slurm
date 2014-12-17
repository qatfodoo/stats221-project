#!/bin/bash
#SBATCH -J ebola-mle # name for job array
#SBATCH -o ./dump/mle.out #Standard output
#SBATCH -e ./dump/mle.err #Standard error
#SBATCH -p general #Partition
#SBATCH -t 01:00:00 #Running time of 1 hours.
#SBATCH --mem-per-cpu 10000 #Memory request
#SBATCH -n 1 #Number of cores
#SBATCH -N 1 #All cores on one machine
#SBATCH --mail-user=kuatefodouop@fas.harvard.edu
#SBATCH --mail-type=ALL    # When to send mail

# argument = job id

Rscript mle_task.R $SLURM_ARRAY_TASK_ID
