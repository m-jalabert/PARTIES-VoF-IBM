#!/bin/bash -l
#SBATCH --job-name=Yang_checkout
#SBATCH --output=checkout_job.out
#SBATCH --error=checkout_job.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=01:00:00
#SBATCH --mail-user=mjalabert@ucsb.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

module load gnu14/14.2.0 hwloc/2.12.0 ucx/1.18.0 libfabric/1.18.0 openmpi5/5.0.7 hdf5/1.14.6

root="$PWD"

run_case() {
    name="$1"
    ranks="$2"
    input="$3"
    mkdir -p "$name"
    cp parties p_fixed.inp p_mobile.inp stop.inp "$name"/
    cp "$input" "$name/parties.inp"
    cd "$name"
    echo "Starting $name on $ranks rank(s)"
    mpirun -np "$ranks" ./parties > run.log 2>&1
    echo "Completed $name"
    cd "$root"
}

run_case run_nz4_1rank 1 parties_nz4.inp
run_case run_nz4_4rank 4 parties_nz4.inp
run_case run_nz8_4rank 4 parties_nz8.inp
run_case run_static_control 4 parties_control.inp

echo "All Yang checkout runs completed"
