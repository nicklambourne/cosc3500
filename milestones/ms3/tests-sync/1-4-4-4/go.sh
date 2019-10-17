#!/bin/bash
#SBATCH --partition=coursework
#SBATCH --job-name=1-4-4-4
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=4
#SBATCH --time=1:00:00

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export TIMEFORMAT="%E sec"

module load mpi/openmpi-x86_64
time mpirun -n ${SLURM_NPROCS} /home/Student/s4261833/cosc3500/milestones/ms3/bin/lcs-hybrid /home/Student/s4261833/cosc3500/milestones/ms3/test/xlong_in.txt /home/Student/s4261833/cosc3500/milestones/ms3/tests/1-4-4-4/1-4-4-4.out

DATE=$(date +"%Y%m%d%H%M")
echo "time finished "$DATE