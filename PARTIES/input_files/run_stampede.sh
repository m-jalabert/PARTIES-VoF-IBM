#!/bin/sh
#SBATCH -A TG-CTS150053
#SBATCH -J g10d_St05wet
#SBATCH -o log.out
#SBATCH -e log.err
#SBATCH -n 32
#SBATCH -t 3:00:00
#SBATCH -p normal
##SBATCH -t 1:00:00
##SBATCH -p development
#SBATCH --mail-user=ebiegert@engineering.ucsb.edu
#SBATCH --mail-user=begin
#SBATCH --mail-user=end

NCPU=$SLURM_NTASKS
NNODES=$SLURM_NNODES

BASDIR=$HOME
CURDIR=$SLURM_SUBMIT_DIR
#RUNDIR=`echo $CURDIR | sed 's|/global/homes/e|/scratch/scratchdirs|'`
RUNDIR=`echo $CURDIR | sed "s|$BASDIR|$WORK|"`
CMD="ibrun"

echo "--> Running on nodes " $SLURM_NODELIST
echo "--> Number of available cpus " $NCPU
echo "--> Number of available nodes " $NNODES
echo "--> Launch command is " $CMD

cd $CURDIR
cp -p parties *.inp $RUNDIR

cd $RUNDIR
$CMD ./parties
