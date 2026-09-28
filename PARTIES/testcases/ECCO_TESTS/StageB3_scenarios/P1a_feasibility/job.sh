#!/bin/bash -l
#SBATCH --job-name=ECCO_P1_feas
#SBATCH --account=phy250167
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --time=03:00:00
#SBATCH --output=/anvil/scratch/x-mjalabert/ECCO_StageB3/P1_Attempts/job_%j.out
#SBATCH --error=/anvil/scratch/x-mjalabert/ECCO_StageB3/P1_Attempts/job_%j.err
set -euo pipefail
module purge
module load gcc/11.2.0 openmpi/4.1.6 hdf5/1.10.7 fftw/3.3.8
build_dir=/anvil/scratch/x-mjalabert/ECCO_StageB3/P1_Attempts/B3_P1_build_20260925
run_dir=/anvil/scratch/x-mjalabert/ECCO_StageB3/B3_P1_20260925_${SLURM_JOB_ID}
mkdir "$run_dir"
cp "$SLURM_SUBMIT_DIR"/*.inp "$SLURM_SUBMIT_DIR"/*.py "$SLURM_SUBMIT_DIR"/*.json "$SLURM_SUBMIT_DIR"/*.md "$SLURM_SUBMIT_DIR"/Boundary.scenario.h "$SLURM_SUBMIT_DIR"/job.sh "$run_dir/"
cp "$build_dir/parties" "$build_dir/manifest.json" "$build_dir/build.log" "$run_dir/"
cp /home/x-mjalabert/PARTIES/PARTIES/testcases/ECCO_TESTS/StageB2_closure/B2_closure/analyze_profiles.py "$run_dir/"
cd "$run_dir"
sha256sum parties *.inp Boundary.scenario.h case.json > run.sha256
module list > modules.txt 2>&1
scontrol show job "$SLURM_JOB_ID" > slurm_start.txt
# Stop through the solver's normal checkpoint path before the hard wall limit.
# Do not resume automatically if melt-out has not occurred.
( sleep 10200; printf '1\n' > stop.inp ) &
watchdog=$!
trap 'kill "$watchdog" 2>/dev/null || true' EXIT
srun ./parties > run.log 2>&1
python analyze_profiles.py . --output budgets > budget_stdout.json
