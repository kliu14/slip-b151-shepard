#!/bin/bash

#SBATCH -J b151-m2
#SBATCH -N 6
#SBATCH -n 192
#SBATCH -p cpu
#SBATCH -o log-%j.out
#SBATCH -e log-%j.err
#SBATCH -t 1-00:00:00
#SBATCH --qos=normal
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kimliu@stanford.edu


WORKDIR='/home/mani/kimliu/slip-shep/'

echo The master node of this job is `hostname`
echo This job runs on the following nodes:
echo `scontrol show hostname $SLURM_JOB_NODELIST`
echo "Starting at `date`"
echo "Running on hosts: $SLURM_NODELIST"
echo "Running on $SLURM_NNODES nodes."
echo "Running on $SLURM_NPROCS processors."
echo "Current working directory is `echo $WORKDIR`"

cd $WORKDIR
mpiexec bin/channel
