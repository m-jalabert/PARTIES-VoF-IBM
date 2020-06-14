#!/bin/sh

#PBS -m be
#PBS -M ebiegert@engineering.ucsb.edu
#PBS -N gondret2_20D
#PBS -o log.out
#PBS -e log.err
#PBS -l mppwidth=120
##PBS -l walltime=4:00:00
##PBS -q regular
#PBS -l walltime=0:30:00
#PBS -q debug

NCPU=`wc -l < $PBS_NODEFILE`
NNODES=`uniq $PBS_NODEFILE | wc -l`

# Establish run directory
BASDIR=$PWD
CURDIR=$PBS_O_WORKDIR
RUNDIR=`echo $CURDIR | sed "s|$BASDIR|$SCRATCH|"`
CMD="aprun -n $NCPU "

echo "--> Running on nodes " `uniq $PBS_NODEFILE`
echo "--> Number of available cpus " $NCPU
echo "--> Number of available nodes " $NNODES
echo "--> Launch command is " $CMD

# Copy input files to run directory
cd $CURDIR
cp -p gvg *.inp $RUNDIR

# Run in run directory
cd $RUNDIR
$CMD ./parties
