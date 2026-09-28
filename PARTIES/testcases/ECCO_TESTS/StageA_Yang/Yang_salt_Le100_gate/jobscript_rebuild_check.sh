#!/bin/bash -l
#SBATCH --job-name=rebuildChk
#SBATCH --output=rebuild_chk.out
#SBATCH --error=rebuild_chk.err
#SBATCH --account=phy250167
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=00:20:00
# Configuration control: the tree is not byte-reproducible (324 bytes differ
# across rebuilds of identical source).  Establish FUNCTIONAL equivalence
# instead, against the recorded legacy gate value lambda = 0.21844955456451656.
set -euo pipefail
module purge
module load gcc/11.2.0 openmpi/4.1.6 hdf5/1.10.7 fftw/3.3.8
C="${SLURM_SUBMIT_DIR}"
d="$C/run_rebuild_chk"; mkdir -p "$d"
cp "$C/parties_legacy.inp" "$d/parties.inp"
cp "$C/p_mobile.inp" "$C/p_fixed.inp" "$C/stop.inp" "$d/"
cp /home/x-mjalabert/PARTIES/PARTIES/parties "$d/parties"
cd "$d" && srun -n 4 ./parties > run.log 2>&1
echo done
