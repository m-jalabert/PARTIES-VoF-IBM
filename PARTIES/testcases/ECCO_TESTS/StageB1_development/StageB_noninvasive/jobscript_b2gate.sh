#!/bin/bash -l
#SBATCH --job-name=B2gate
#SBATCH --output=b2gate.out
#SBATCH --error=b2gate.err
#SBATCH --account=phy250167
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=00:40:00
# B2 default-off regression: the interface-release code must be a no-op when
# F_release <= 0 (its default).  Same deck, same reference, new binary.
set -euo pipefail
module purge
module load gcc/11.2.0 openmpi/4.1.6 hdf5/1.10.7 fftw/3.3.8
C="${SLURM_SUBMIT_DIR}"
d="$C/gate_b2off"; mkdir -p "$d"
cp "$C/parties_legacy.inp" "$d/parties.inp"
cp "$C/stop.inp" "$C/p_fixed.inp" "$C/p_mobile.inp" "$d/"
cp /home/x-mjalabert/PARTIES/PARTIES/parties.stageB_release "$d/parties"
cd "$d" && srun ./parties > run.log 2>&1
echo "b2 gate done"
