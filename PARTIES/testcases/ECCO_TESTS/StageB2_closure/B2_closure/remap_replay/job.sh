#!/bin/bash -l
#SBATCH --job-name=B2remap_replay
#SBATCH --account=phy250167
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=00:30:00
#SBATCH --output=job_%j.out
#SBATCH --error=job_%j.err
set -eo pipefail
module purge
module load gcc/11.2.0 openmpi/4.1.6 hdf5/1.10.7 fftw/3.3.8
run_dir=/anvil/scratch/x-mjalabert/ECCO_StageB/B2closure_remap_replay_${SLURM_JOB_ID}
mkdir -p "$run_dir"
cp "$SLURM_SUBMIT_DIR"/*.inp "$run_dir/"
cp /home/x-mjalabert/PARTIES/PARTIES/testcases/ECCO_TESTS/StageB2_closure/B2_closure/parties.ice_v8 "$run_dir/parties"
sha256sum "$run_dir/parties" > "$run_dir/binary.sha256"
cp /anvil/scratch/x-mjalabert/ECCO_StageB/B2closure_coupled_v6_20477118/Resume_25.h5 "$run_dir/Resume.h5"
cp /anvil/scratch/x-mjalabert/ECCO_StageB/B2closure_coupled_v6_20477118/Data_25.h5 /anvil/scratch/x-mjalabert/ECCO_StageB/B2closure_coupled_v6_20477118/Particle_25.h5 "$run_dir/"
cd "$run_dir"
srun ./parties > run.log 2>&1
