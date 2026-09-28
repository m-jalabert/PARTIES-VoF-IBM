#!/bin/bash -l
#SBATCH --job-name=sliqGate
#SBATCH --output=sliq_job.out
#SBATCH --error=sliq_job.err
#SBATCH --account=phy250167
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=00:30:00

# A.5.5 variant 1 screen.  sliq0 MUST reproduce the recorded legacy
# lambda_fit = 0.21844955456451656 -- that is the regression test proving the
# new binary is unchanged with the flag off.  sliq1 is the variant.
# NOTE: this gate is quiescent with eos_betaT = eos_betaS = 0, so it tests ONLY
# the liquidus half -- but at liquidus_slope = 0.5, i.e. 36x the production
# value 0.014, so it is a sensitive screen for that half.  The EOS half
# (betaS = 1.75) can only be tested in 2-D.
set -euo pipefail
module purge
module load gcc/11.2.0 openmpi/4.1.6 hdf5/1.10.7 fftw/3.3.8
CASE_DIR="${SLURM_SUBMIT_DIR}"
PARTIES_BIN="${PARTIES_BIN:-/home/x-mjalabert/PARTIES/PARTIES/parties}"
for b in sliq0 sliq1; do
    d="${CASE_DIR}/run_${b}"; mkdir -p "$d"
    cp "${CASE_DIR}/parties_${b}.inp" "$d/parties.inp"
    cp "${CASE_DIR}/p_mobile.inp" "${CASE_DIR}/p_fixed.inp" "${CASE_DIR}/stop.inp" "$d/"
    cp "${PARTIES_BIN}" "$d/parties"
    cd "$d"; srun -n 4 ./parties > run.log 2>&1; cd "${CASE_DIR}"
    echo "$b done"
done
