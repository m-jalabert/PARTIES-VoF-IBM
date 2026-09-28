#!/bin/bash -l
#SBATCH --job-name=StageBnoninv
#SBATCH --output=noninv.out
#SBATCH --error=noninv.err
#SBATCH --account=phy250167
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --time=01:00:00
# Does enabling VOF_IBM + LAG_PARTICLE_RESOLVED perturb the validated Stage-A
# physics when no particles are present?  The roadmap's non-invasiveness premise
# has never been tested with these flags ON.
set -euo pipefail
module purge
module load gcc/11.2.0 openmpi/4.1.6 hdf5/1.10.7 fftw/3.3.8
C="${SLURM_SUBMIT_DIR}"; B=/home/x-mjalabert/PARTIES/PARTIES/parties.stageB

d="$C/gate_stageB"; mkdir -p "$d"
cp "$C/parties_legacy.inp" "$d/parties.inp"; cp "$C/p_mobile.inp" "$C/p_fixed.inp" "$C/stop.inp" "$d/"
cp "$B" "$d/parties"; cd "$d"; srun -n 4 ./parties > run.log 2>&1; echo "gate done"; cd "$C"

d="$C/fresh288_stageB"; mkdir -p "$d"
cp "$C/parties_fresh288.inp" "$d/parties.inp"; cp "$C/stop.inp" "$d/"
cp "$B" "$d/parties"; cd "$d"; srun ./parties > run.log 2>&1; echo "fresh288 done"
