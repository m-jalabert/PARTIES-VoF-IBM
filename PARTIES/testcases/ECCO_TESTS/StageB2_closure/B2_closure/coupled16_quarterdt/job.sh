#!/bin/bash -l
#SBATCH --job-name=B2coupled16_quarterdt
#SBATCH --account=phy250167
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --time=06:00:00
#SBATCH --output=job_%j.out
#SBATCH --error=job_%j.err
set -eo pipefail
module purge
module load gcc/11.2.0 openmpi/4.1.6 hdf5/1.10.7 fftw/3.3.8
run_dir=/anvil/scratch/x-mjalabert/ECCO_StageB/B2closure_coupled16_quarterdt_${SLURM_JOB_ID}
mkdir -p "$run_dir"
cp "$SLURM_SUBMIT_DIR"/*.inp "$run_dir/"
cp /home/x-mjalabert/PARTIES/PARTIES/testcases/ECCO_TESTS/StageB2_closure/B2_closure/parties.ice_v9 "$run_dir/parties"
sha256sum "$run_dir/parties" > "$run_dir/binary.sha256"
cd "$run_dir"
srun ./parties > run.log 2>&1
