#!/bin/bash -l
#SBATCH --job-name=B14rfctl
#SBATCH --output=rfctl.out
#SBATCH --error=rfctl.err
#SBATCH --account=phy250167
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=06:00:00
#
# Roadmap B.1.4 -- CONTROL arm.  The StageB_release2d deck, unchanged, run with
# SED_CLIP_AUDIT on and the sediment-clip restore OFF.  This is the roadmap's
# "instrument first, decide on data" step: it measures the liquid destroyed in
# the ice/sediment overlap directly (SED_CLIP_AUDIT lines in run.log) instead
# of inferring it from the meltwater-tracer budget, which conflates it with
# tracer transport error.
#
# Pairs with jobscript_fix.sh, which differs by exactly one compile flag.
set -euo pipefail
set +u; module purge; module load gcc/11.2.0 openmpi/4.1.6 hdf5/1.10.7 fftw/3.3.8; set -u
C="${SLURM_SUBMIT_DIR}"
# Run output goes to SCRATCH, never $HOME (25 GB quota).  See roadmap B.0.1.
d="/anvil/scratch/x-mjalabert/ECCO_StageB/B14_relfn_ctl_${SLURM_JOB_ID}"
mkdir -p "$d"
cp "$C/parties.inp" "$C/stop.inp" "$C/p_fixed.inp" "$C/p_mobile.inp" "$d/"
cp /home/x-mjalabert/PARTIES/PARTIES/parties.stageB_sedctl "$d/parties"
cd "$d" && srun ./parties > run.log 2>&1
echo "B.1.4 control arm done"
