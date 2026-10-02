#!/bin/bash -l
# P1b follow-up variants (RESULTS_20977443.md sec. 6).  Usage, from a staged copy:
#   sbatch --export=ALL,VARIANT=conc4,BUILD=build_v2 --time=12:00:00 job_variant.sh
#   sbatch --export=ALL,VARIANT=twin_lid_strat,BUILD=build_v3 --time=04:00:00 job_variant.sh
# The variant deck and p_mobile.inp come from variants/<VARIANT>/; output goes to a
# job-unique directory.  No automatic continuation.
#SBATCH --job-name=ECCO_P1b_var
#SBATCH --account=phy250167
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --time=12:00:00
#SBATCH --output=/anvil/scratch/x-mjalabert/ECCO_StageB3/P1_Attempts/P1b_shelf/logs/var_%j.out
#SBATCH --error=/anvil/scratch/x-mjalabert/ECCO_StageB3/P1_Attempts/P1b_shelf/logs/var_%j.err
set -euo pipefail
: "${VARIANT:?set VARIANT}" "${BUILD:?set BUILD}"
module purge
module load gcc/11.2.0 openmpi/4.1.6 hdf5/1.10.7 fftw/3.3.8
ROOT=/anvil/scratch/x-mjalabert/ECCO_StageB3/P1_Attempts/P1b_shelf
run_dir=$ROOT/run_${VARIANT}_${SLURM_JOB_ID}
mkdir "$run_dir"
cp "$SLURM_SUBMIT_DIR"/variants/"$VARIANT"/parties.inp "$SLURM_SUBMIT_DIR"/variants/"$VARIANT"/p_mobile.inp \
   "$SLURM_SUBMIT_DIR"/p_fixed.inp "$SLURM_SUBMIT_DIR"/stop.inp "$SLURM_SUBMIT_DIR"/case.json \
   "$SLURM_SUBMIT_DIR"/Boundary.scenario.h "$SLURM_SUBMIT_DIR"/*.py "$SLURM_SUBMIT_DIR"/job_variant.sh "$run_dir/"
cp "$ROOT/$BUILD/parties" "$ROOT/$BUILD/manifest.json" "$run_dir/"
cp /home/x-mjalabert/PARTIES/PARTIES/testcases/ECCO_TESTS/StageB2_closure/B2_closure/analyze_profiles.py "$run_dir/"
cd "$run_dir"
echo "$VARIANT $BUILD" > case_label.txt
sha256sum parties *.inp Boundary.scenario.h case.json > run.sha256
scontrol show job "$SLURM_JOB_ID" > slurm_start.txt
# checkpoint-stop 20 min before the requested wall time
limit=$(squeue -h -j "$SLURM_JOB_ID" -o %L | awk -F'[-:]' '{n=NF; s=$n+60*$(n-1)+3600*$(n-2); if(n==4) s+=86400*$1; print s}')
( sleep $(( limit > 1500 ? limit - 1200 : limit / 2 )); printf '1\n' > stop.inp ) &
watchdog=$!
trap 'kill "$watchdog" 2>/dev/null || true' EXIT
srun ./parties > run.log 2>&1
python analyze_profiles.py . --output budgets > budget_stdout.json || true
