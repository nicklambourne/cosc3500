#!/bin/bash
#SBATCH --partition=coursework
#SBATCH --job-name={{name}}
#SBATCH --nodes={{nnodes}}
#SBATCH --ntasks={{ntasks}}
#SBATCH --ntasks-per-node={{ntasks_per_node}}
#SBATCH --cpus-per-task={{cpus_per_task}}

# export SLURM_NNODES=1
# export SLURM_NTASKS=4
# export SLURM_TASKS_PER_NODE=4
# export SLURM_CPUS_PER_TASK=4
# export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

DATE=$(date +"%Y%m%d%H%M")
echo "time started "$DATE
echo "This is job '$SLURM_JOB_NAME' (id: $SLURM_JOB_ID) running on the following nodes:"
echo $SLURM_NODELIST
echo "running with OMP_NUM_THREADS= $OMP_NUM_THREADS "
echo "running with SLURM_TASKS_PER_NODE= $SLURM_TASKS_PER_NODE "
echo "running with SLURM_NPROCS= $SLURM_NPROCS "
echo "Now we start the show:"
export TIMEFORMAT="%E sec"

module load mpi/openmpi-x86_64
time mpirun -n ${SLURM_NPROCS} ./bin/lcs-hybrid ./test/xlong_in.txt ./test/xlong_out.out

DATE=$(date +"%Y%m%d%H%M")
echo "time finished "$DATE
# echo "we just ran with the following SLURM environment variables"
# echo $OMP_NUM_THREADS
# env | grep SLURM