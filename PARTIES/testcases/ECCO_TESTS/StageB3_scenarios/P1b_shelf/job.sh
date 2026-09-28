#!/bin/bash -l
# P1b production: released fine-sand grain beneath an Antarctic ice-shelf base.
#   sbatch job.sh                      -> P1b (grain)
#   sbatch --export=ALL,CASE=ctrl job.sh -> matched ice-only control (no grain)
# Submit from a staged copy of this folder (never from home); output goes to a
# job-unique run directory.  No automatic continuation.
#SBATCH --job-name=ECCO_P1b_shelf
#SBATCH --account=phy250167
#SBATCH --partition=wholenode
#SBATCH --nodes=2
#SBATCH --ntasks=256
#SBATCH --time=12:00:00
#SBATCH --output=/anvil/scratch/x-mjalabert/ECCO_StageB3/P1_Attempts/P1b_shelf/logs/P1b_%j.out
#SBATCH --error=/anvil/scratch/x-mjalabert/ECCO_StageB3/P1_Attempts/P1b_shelf/logs/P1b_%j.err
set -euo pipefail
CASE=${CASE:-prod}
module purge
module load gcc/11.2.0 openmpi/4.1.6 hdf5/1.10.7 fftw/3.3.8
ROOT=/anvil/scratch/x-mjalabert/ECCO_StageB3/P1_Attempts/P1b_shelf
run_dir=$ROOT/run_${CASE}_${SLURM_JOB_ID}
mkdir "$run_dir"
cp "$SLURM_SUBMIT_DIR"/parties.inp "$SLURM_SUBMIT_DIR"/p_mobile.inp "$SLURM_SUBMIT_DIR"/p_fixed.inp \
   "$SLURM_SUBMIT_DIR"/stop.inp "$SLURM_SUBMIT_DIR"/case.json "$SLURM_SUBMIT_DIR"/Boundary.scenario.h \
   "$SLURM_SUBMIT_DIR"/*.py "$SLURM_SUBMIT_DIR"/job.sh "$run_dir/"
if [ "$CASE" = ctrl ]; then
    cp "$SLURM_SUBMIT_DIR"/ctrl/parties.inp "$SLURM_SUBMIT_DIR"/ctrl/p_mobile.inp "$run_dir/"
fi
cp "$ROOT/build/parties" "$ROOT/build/manifest.json" "$ROOT/build/build.log" "$run_dir/"
cp /home/x-mjalabert/PARTIES/PARTIES/testcases/ECCO_TESTS/StageB2_closure/B2_closure/analyze_profiles.py "$run_dir/"
cd "$run_dir"
echo "$CASE" > case_label.txt
sha256sum parties *.inp Boundary.scenario.h case.json > run.sha256
module list > modules.txt 2>&1
scontrol show job "$SLURM_JOB_ID" > slurm_start.txt
# Checkpoint-stop through the solver's normal path 30 min before the wall limit.
( sleep 41400; printf '1\n' > stop.inp ) &
watchdog=$!
trap 'kill "$watchdog" 2>/dev/null || true' EXIT
srun ./parties > run.log 2>&1
python analyze_profiles.py . --output budgets > budget_stdout.json || true
