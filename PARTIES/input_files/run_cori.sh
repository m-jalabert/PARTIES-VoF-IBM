#!/bin/bash -l

#SBATCH --mail-type=ALL
##SBATCH --mail-users=ebiegert@engineering.ucsb.edu
#SBATCH -J gondret10d_sense_00
#SBATCH -o log.out
#SBATCH -e log.err
#SBATCH -N 16
##SBATCH -p regular
##SBATCH -t 06:00:00
#SBATCH -p debug
#SBATCH -t 00:30:00

NCPU=$((32*$SLURM_JOB_NUM_NODES))

# Establish run directory
BASDIR=$HOME
CURDIR=$PWD
RUNDIR=`echo $CURDIR | sed "s|$BASDIR|$SCRATCH|"`
CMD="srun -n $NCPU "

echo "--> Running on nodes " $SLURM_NODELIST
echo "--> Number of available cpus " $NCPU
echo "--> Number of available nodes " $SLURM_JOB_NUM_NODES
echo "--> Launch command is " $CMD

# Copy input files to run directory
cd $CURDIR
cp -p parties *.inp $RUNDIR

# Run in run directory
cd $RUNDIR
$CMD ./parties

