#!/bin/bash -l
#SBATCH --job-name=yang_salt_retry
#SBATCH --output=gate_retry_job.out
#SBATCH --error=gate_retry_job.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=00:10:00

set -euo pipefail

module purge
module load gcc/11.2.0 openmpi/4.1.6 hdf5/1.10.7 fftw/3.3.8

CASE_DIR="${SLURM_SUBMIT_DIR}"
PARTIES_BIN="${PARTIES_BIN:-/home/x-mjalabert/PARTIES/PARTIES/parties}"
run_dir="${CASE_DIR}/run_yang"

mkdir -p "${run_dir}"
cp "${CASE_DIR}/parties_yang.inp" "${run_dir}/parties.inp"
cp "${CASE_DIR}/p_mobile.inp" "${CASE_DIR}/p_fixed.inp" "${CASE_DIR}/stop.inp" "${run_dir}/"
cp "${PARTIES_BIN}" "${run_dir}/parties"

cd "${run_dir}"
mpirun -np "${SLURM_NTASKS}" ./parties > run.log 2>&1

cd "${CASE_DIR}"
python3 analyze_gate.py --legacy run_legacy --yang run_yang
