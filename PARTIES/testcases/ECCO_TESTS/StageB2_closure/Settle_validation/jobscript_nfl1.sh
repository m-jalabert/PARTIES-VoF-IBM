#!/bin/bash -l
#SBATCH --job-name=SETTLE1
#SBATCH --output=settle1.out
#SBATCH --error=settle1.err
#SBATCH --account=phy250167
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --time=03:00:00
#
# Roadmap B.2 -- reduced 3-D shakeout.  First time the FULL coupled path runs in
# 3-D: melt + motion lock + interface-triggered release + settling + three
# scalars + bottom-wall collision.  Code shake-out, NOT science (d/dx = 16,
# Sc reduced).
#
# Binary parties.stageB_sed3d:
#   * x,z periodic + y no-slip (PERIODIC_NOSLIP_BOX), TWOD_CARTESIAN off
#   * VOF_DIFFUSE_SEDIMENT_CLIP_RESTORE on  (B.1.4 fix -- the leak grows with
#     sediment surface area, so 3-D is where it matters most)
#   * VOF_DIFFUSE_SEDIMENT_CLIP_AUDIT 200   (measures it while we are here)
#   * LUBRICATION_NORMAL on, so the collision path exercised here is the one
#     B.3 production will actually use
#
# Deck change vs the Jul-30 draft: F_release 0.9 -> 0.709, the B.1.3-fn
# calibration ("released the moment it is mechanically free"); 0.9 fires 0.67
# diameters late.
#
# What to check (B.2): (1) particle exactly stationary while locked;
# (2) impulse-free release, peak |dU/dt| < 2x free-settling; (3) meltwater
# tracer injected only at the interface and globally conserved; (4) wall
# collision fires without blow-up; (5) heat/salt/tracer budgets close.
set -euo pipefail
set +u; module purge; module load gcc/11.2.0 openmpi/4.1.6 hdf5/1.10.7 fftw/3.3.8; set -u
C="${SLURM_SUBMIT_DIR}"
# Run output goes to SCRATCH, never $HOME (25 GB quota).  See roadmap B.0.1.
d="/anvil/scratch/x-mjalabert/ECCO_StageB/Settle_nfl1_${SLURM_JOB_ID}"
mkdir -p "$d"
cp "$C/parties.inp" "$C/stop.inp" "$C/p_fixed.inp" "$C/p_mobile.inp" "$d/"
cp /home/x-mjalabert/PARTIES/PARTIES/parties.stageB_sed3d "$d/parties"
cd "$d" && srun ./parties > run.log 2>&1
echo "B.2 3-D shakeout done"
