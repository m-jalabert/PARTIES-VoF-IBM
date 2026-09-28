#!/bin/bash -l
#SBATCH --job-name=B2restart_guards
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
run_dir=/anvil/scratch/x-mjalabert/ECCO_StageB/B2closure_restart_guards_${SLURM_JOB_ID}
mkdir -p "$run_dir"
cp "$SLURM_SUBMIT_DIR"/*.inp "$run_dir/"
cp /home/x-mjalabert/PARTIES/PARTIES/testcases/ECCO_TESTS/StageB2_closure/B2_closure/parties.ice_v9 "$run_dir/parties"
sha256sum "$run_dir/parties" > "$run_dir/binary.sha256"
for test_case in mismatch missing_data; do
  mkdir -p "$run_dir/$test_case"
  cp "$run_dir"/*.inp "$run_dir/$test_case/"
  if [[ "$test_case" == mismatch ]]; then
    cp /anvil/scratch/x-mjalabert/ECCO_StageB/B2closure_profile_seed_20477117/{Resume.h5,Data_2.h5,Particle_2.h5} "$run_dir/$test_case/"
    expected=93
    message="Restart phase transport mismatch"
  else
    cp /anvil/scratch/x-mjalabert/ECCO_StageB/B2closure_seed_v9_20477869/{Resume.h5,Particle_2.h5} "$run_dir/$test_case/"
    expected=92
    message="Cannot open restart file"
  fi
  cd "$run_dir/$test_case"
  set +e
  srun ../parties > run.log 2>&1
  actual=$?
  set -e
  if [[ "$actual" == 0 ]] || ! rg -q "$message" run.log || ! rg -q "errorcode $expected\.|PMI_Abort\($expected," run.log; then
    echo "FAIL $test_case exit=$actual expected=$expected"
    exit 1
  fi
  echo "PASS $test_case launcher_exit=$actual MPI_Abort=$expected"
done
