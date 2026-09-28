#!/bin/bash -l
#SBATCH --job-name=sliqBisect
#SBATCH --output=bisect.out
#SBATCH --error=bisect.err
#SBATCH --account=phy250167
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=00:30:00
set -uo pipefail
module purge
module load gcc/11.2.0 openmpi/4.1.6 hdf5/1.10.7 fftw/3.3.8
C="${SLURM_SUBMIT_DIR}"; P=/home/x-mjalabert/PARTIES/PARTIES
# 1 = control: reference binary, original input
# 2 = isolates the code change: new binary, original input (no new key)
# 3 = flag off; 4 = flag on
run () {  # name binary input
  d="$C/bis_$1"; mkdir -p "$d"
  cp "$C/$3" "$d/parties.inp"
  cp "$C/p_mobile.inp" "$C/p_fixed.inp" "$C/stop.inp" "$d/"
  cp "$2" "$d/parties"
  cd "$d" && srun -n 4 ./parties > run.log 2>&1
  echo "$1 exit=$?"
  cd "$C"
}
run ref_orig  "$P/parties.ref_6a8fa6eb" parties_legacy.inp
run new_orig  "$P/parties"              parties_legacy.inp
run new_off   "$P/parties"              parties_sliq0.inp
run new_on    "$P/parties"              parties_sliq1.inp
