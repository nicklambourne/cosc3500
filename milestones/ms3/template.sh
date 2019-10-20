#!/bin/bash
#SBATCH --partition=coursework
#SBATCH --job-name={{name}}
#SBATCH --nodes={{nnodes}}
#SBATCH --ntasks={{ntasks}}
#SBATCH --ntasks-per-node={{ntasks_per_node}}
#SBATCH --cpus-per-task={{cpus_per_task}}
#SBATCH --time=1:00:00

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export TIMEFORMAT="%E sec"

echo "{{name}} - {{input}}"

module load mpi/openmpi-x86_64
time mpirun -n ${SLURM_NPROCS} {{binary}} {{input}} {{output}}

DATE=$(date +"%Y%m%d%H%M")
echo "time finished "$DATE