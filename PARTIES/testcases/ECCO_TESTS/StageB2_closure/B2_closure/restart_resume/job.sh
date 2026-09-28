#!/bin/bash -l
#SBATCH --job-name=B2restart_resume
#SBATCH --account=phy250167
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=00:15:00
#SBATCH --output=job_%j.out
#SBATCH --error=job_%j.err
set -eo pipefail
module purge
module load gcc/11.2.0 openmpi/4.1.6 hdf5/1.10.7 fftw/3.3.8
run_dir=/anvil/scratch/x-mjalabert/ECCO_StageB/B2closure_restart_resume_${SLURM_JOB_ID}
mkdir -p "$run_dir"
cp "$SLURM_SUBMIT_DIR"/*.inp "$run_dir/"
cp /home/x-mjalabert/PARTIES/PARTIES/testcases/ECCO_TESTS/StageB2_closure/B2_closure/parties.ice_v7 "$run_dir/parties"
sha256sum "$run_dir/parties" > "$run_dir/binary.sha256"
cp /anvil/scratch/x-mjalabert/ECCO_StageB/B2closure_profile_seed_20477117/Resume.h5 /anvil/scratch/x-mjalabert/ECCO_StageB/B2closure_profile_seed_20477117/Data_2.h5 /anvil/scratch/x-mjalabert/ECCO_StageB/B2closure_profile_seed_20477117/Particle_2.h5 "$run_dir/"
cd "$run_dir"
srun ./parties > run.log 2>&1
