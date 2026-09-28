#!/bin/bash -l
#SBATCH --job-name=SB2_3dR
#SBATCH --output=sb3d_resume.out
#SBATCH --error=sb3d_resume.err
#SBATCH --account=phy250167
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --time=24:00:00
#
# Roadmap B.2 -- RESUME of the 3-D shakeout.
#
# Job 20218407 runs at ~9.0 t.u./h, so its 24 h wall reaches t ~ 216 of
# time_max = 250.  By the one-phase Stefan law the front needs t ~ 200 to clear
# the grain (X = 0.0708*sqrt(t) at Pe_T = 100, 1.0 length unit of ice above the
# grain), so the release lands right at the edge of that window -- too tight to
# rely on.  ~28 h are needed for the full run.
#
# Usage:  set PREV to the timed-out run directory, then submit.  The deck is
# copied from it (not from $SLURM_SUBMIT_DIR) so the resumed run uses exactly
# the inputs the original ran with, with resume flipped on.
set -euo pipefail
set +u; module purge; module load gcc/11.2.0 openmpi/4.1.6 hdf5/1.10.7 fftw/3.3.8; set -u

PREV="${PREV:-/anvil/scratch/x-mjalabert/ECCO_StageB/StageB_release3d_20240184}"
TIME_MAX="${TIME_MAX:-450.0}"
d="/anvil/scratch/x-mjalabert/ECCO_StageB/StageB_release3d_resume_${SLURM_JOB_ID}"
mkdir -p "$d"

cp "$PREV"/{parties.inp,stop.inp,p_fixed.inp,p_mobile.inp} "$d/"
# Carry the restart state.  Resume.h5 is a symlink to the newest Resume_N.h5,
# so dereference it (-L) rather than copying a dangling link.
cp -L "$PREV/Resume.h5" "$d/Resume.h5"
cp -L "$PREV"/Particle_*.h5 "$d/" 2>/dev/null || true
# Resume.c has TWO read paths: one opens Resume.h5 (line 67) and a second opens
# Data_<noutput>.h5 (line 382) for the Eulerian fields.  Staging only Resume.h5
# gives H5Fopen a bad file_id, which is NOT checked -- the run dies with
# 'Open dataset "/u" from <garbage> failed' and a SIGSEGV (seen in job 20253869).
# Carry the last two Data snapshots so whichever index noutput holds is present.
for n in $(ls "$PREV"/Data_*.h5 | sed 's/.*Data_//; s/\.h5//' | sort -n | tail -2); do
    cp -L "$PREV/Data_${n}.h5" "$d/"
done
sed -i 's/^resume = 0$/resume = 1/' "$d/parties.inp"
# Extend the run: at t=250 phi_liq = 0.6436 against the 0.709 threshold, rising a
# steady ~0.038 per 25 t.u. (the warm bottom wall made melting sustained), so the
# release lands near t = 293.  The extra headroom is for the settle + wall collision
# that B.2 checks (2) and (4) actually need.
sed -i "s/^time_max = .*/time_max = ${TIME_MAX}/" "$d/parties.inp"

cp /home/x-mjalabert/PARTIES/PARTIES/parties.stageB_sed3d "$d/parties"
cd "$d" && srun ./parties > run.log 2>&1
echo "B.2 3-D shakeout resume done"
