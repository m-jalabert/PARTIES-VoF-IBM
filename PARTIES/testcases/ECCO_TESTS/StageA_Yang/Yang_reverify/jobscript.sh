#!/bin/bash -l
#SBATCH --job-name=YangReVer
#SBATCH --output=yrv.out
#SBATCH --error=yrv.err
#SBATCH --account=phy250167
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=02:00:00
#
# Does the B.1.4 / B.2 work change a CONVECTIVE Stage-A result?
#
# The Le=100 no-impact gate has max speed EXACTLY 0, so it cannot exercise the
# velocity solver or ICE_PENALIZATION at all -- and ICE_PENALIZATION *is* defined
# in the Stage-A/Yang/RB flag set, while the 2026-08-31 work edited msolve_cg.c,
# msolve_direct.c and Velocity.c.  The logical argument that this is inert
# (C_S == 0 => cs_face == 0.0 exactly, and x - 0.0 == x in IEEE) is sound but
# untested in a flowing case.  This run tests it.
#
# Yang freshwater N=288 to t=10, against the recorded reference at
# /anvil/scratch/x-mjalabert/YangFresh_07272026/N288.  Standard is B.0 part 2:
# "indistinguishable from a rebuild", not bit-identity -- this case is chaotic
# (Lyapunov rate 0.652/t.u.), so compare ice volume at matching output times.
set -euo pipefail
set +u; module purge; module load gcc/11.2.0 openmpi/4.1.6 hdf5/1.10.7 fftw/3.3.8; set -u
C="${SLURM_SUBMIT_DIR}"
d="/anvil/scratch/x-mjalabert/ECCO_StageB/Yang_reverify_${SLURM_JOB_ID}"
mkdir -p "$d"
cp "$C/parties.inp" "$C/stop.inp" "$d/"
cp /home/x-mjalabert/PARTIES/PARTIES/parties.stageA_b14 "$d/parties"
cd "$d" && srun ./parties > run.log 2>&1
echo "Yang convective re-verification done"
