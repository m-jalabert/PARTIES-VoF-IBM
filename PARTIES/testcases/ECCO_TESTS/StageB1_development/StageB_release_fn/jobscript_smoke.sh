#!/bin/bash -l
#SBATCH --job-name=B13smoke
#SBATCH --output=b13smoke.out
#SBATCH --error=b13smoke.err
#SBATCH --account=phy250167
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=00:30:00
set -euo pipefail
module purge
module load gcc/11.2.0 openmpi/4.1.6 hdf5/1.10.7 fftw/3.3.8
C="${SLURM_SUBMIT_DIR}"
d="$C/smoke"; rm -rf "$d"; mkdir -p "$d"
cp "$C/parties_smoke.inp" "$d/parties.inp"
cp "$C/stop.inp" "$C/p_fixed.inp" "$C/p_mobile.inp" "$d/"
cp /home/x-mjalabert/PARTIES/PARTIES/parties.stageB_rel2 "$d/parties"
cd "$d" && srun ./parties > run.log 2>&1
echo "smoke done"
