#!/bin/bash -l
# P1b pre-flight (debug queue, < 2 h): four 1-D kernel-gate variants plus a
# full-size smoke run of the production deck, all with the isolated P1b build.
# Submit from a staged copy of this folder; output goes to a job-unique dir.
#SBATCH --job-name=P1b_gate_smoke
#SBATCH --account=phy250167
#SBATCH --partition=debug
#SBATCH --nodes=2
#SBATCH --ntasks=256
#SBATCH --time=01:50:00
#SBATCH --output=/anvil/scratch/x-mjalabert/ECCO_StageB3/P1_Attempts/P1b_shelf/logs/gate_smoke_%j.out
#SBATCH --error=/anvil/scratch/x-mjalabert/ECCO_StageB3/P1_Attempts/P1b_shelf/logs/gate_smoke_%j.err
set -uo pipefail
module purge
module load gcc/11.2.0 openmpi/4.1.6 hdf5/1.10.7 fftw/3.3.8
ROOT=/anvil/scratch/x-mjalabert/ECCO_StageB3/P1_Attempts/P1b_shelf
BIN=$ROOT/build/parties
SUB=$SLURM_SUBMIT_DIR
OUT=$ROOT/gate_smoke_${SLURM_JOB_ID}
mkdir "$OUT"
cp -r "$SUB"/gate "$SUB"/case.json "$SUB"/*.py "$SUB"/job_gate_smoke.sh "$OUT/"
cp "$BIN" "$ROOT/build/manifest.json" "$OUT/"
sha256sum "$BIN" > "$OUT/binary.sha256"
scontrol show job "$SLURM_JOB_ID" > "$OUT/slurm_start.txt"

for v in G1_sliq1_n24 G0_sliq0_n24 G1_sliq1_n24_halfdt G1_sliq1_n32; do
    d=$OUT/gate/$v
    printf '0\n' > "$d/p_mobile.inp"; printf '0\n' > "$d/p_fixed.inp"; printf '0\n' > "$d/stop.inp"
    cp "$BIN" "$d/parties"
    ( cd "$d" && if srun --exact -n 8 ./parties > run.log 2>&1; then echo "$v ok"; else echo "$v FAILED"; fi >> "$OUT/status.txt" ) &
done

d=$OUT/smoke
mkdir "$d"
cp "$SUB"/smoke/parties.inp "$SUB"/p_mobile.inp "$SUB"/p_fixed.inp "$SUB"/stop.inp "$d/"
cp "$BIN" "$d/parties"
( cd "$d" && if srun --exact -n 224 ./parties > run.log 2>&1; then echo "smoke ok"; else echo "smoke FAILED"; fi >> "$OUT/status.txt" ) &
wait
echo "all steps finished" >> "$OUT/status.txt"
