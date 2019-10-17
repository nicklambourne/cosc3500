#!/bin/bash
#SBATCH --partition=coursework
#SBATCH --job-name=2-2-1-4
#SBATCH --nodes=2
#SBATCH --ntasks=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=1:00:00

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export TIMEFORMAT="%E sec"

module load mpi/openmpi-x86_64
time mpirun -n ${SLURM_NPROCS} /home/Student/s4261833/cosc3500/milestones/ms3/bin/lcs-hybrid /home/Student/s4261833/cosc3500/milestones/ms3/test/xlong_in.txt /home/Student/s4261833/cosc3500/milestones/ms3/tests/2-2-1-4/2-2-1-4.out

DATE=$(date +"%Y%m%d%H%M")
echo "time finished "$DATE