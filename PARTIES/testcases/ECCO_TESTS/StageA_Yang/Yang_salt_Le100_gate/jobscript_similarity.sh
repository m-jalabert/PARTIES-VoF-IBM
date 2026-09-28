#!/bin/bash -l
#SBATCH --job-name=yang_similarity
#SBATCH --output=similarity_job.out
#SBATCH --error=similarity_job.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=00:10:00

set -euo pipefail

module purge
module load gcc/11.2.0 openmpi/4.1.6 hdf5/1.10.7 fftw/3.3.8

CASE_DIR="${SLURM_SUBMIT_DIR}"
PARTIES_BIN="${PARTIES_BIN:-/home/x-mjalabert/PARTIES/PARTIES/parties}"

cd "${CASE_DIR}"
python3 make_similarity_seed.py --source run_legacy --output similarity_seed

for branch in legacy yang; do
    run_dir="${CASE_DIR}/run_similarity_${branch}"
    mkdir -p "${run_dir}"
    cp "${CASE_DIR}/parties_similarity_${branch}.inp" "${run_dir}/parties.inp"
    cp "${CASE_DIR}/p_mobile.inp" "${CASE_DIR}/p_fixed.inp" "${CASE_DIR}/stop.inp" "${run_dir}/"
    cp "${CASE_DIR}/similarity_seed/Data_5.h5" "${CASE_DIR}/similarity_seed/Resume.h5" "${run_dir}/"
    cp "${PARTIES_BIN}" "${run_dir}/parties"

    cd "${run_dir}"
    mpirun -np "${SLURM_NTASKS}" ./parties > run.log 2>&1
done

cd "${CASE_DIR}"
python3 analyze_gate.py \
    --legacy run_similarity_legacy \
    --yang run_similarity_yang \
    --output-prefix similarity_gate
