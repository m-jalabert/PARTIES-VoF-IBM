#!/bin/bash -l
# P1b pre-flight round 2 (debug queue):
#  * gate round 2: Yang salt operator (G2_yang_n24, G2_yang_n32), production binary build_v2;
#  * Poisson-stall diagnostics D1-D7 (small box, first ~15 steps), DIAGNOSTIC binary
#    build_diag (= build_v2 + divergence-location print), never used for production.
#SBATCH --job-name=P1b_preflight2
#SBATCH --account=phy250167
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --time=01:00:00
#SBATCH --output=/anvil/scratch/x-mjalabert/ECCO_StageB3/P1_Attempts/P1b_shelf/logs/preflight2_%j.out
#SBATCH --error=/anvil/scratch/x-mjalabert/ECCO_StageB3/P1_Attempts/P1b_shelf/logs/preflight2_%j.err
set -uo pipefail
module purge
module load gcc/11.2.0 openmpi/4.1.6 hdf5/1.10.7 fftw/3.3.8
ROOT=/anvil/scratch/x-mjalabert/ECCO_StageB3/P1_Attempts/P1b_shelf
SUB=$SLURM_SUBMIT_DIR
OUT=$ROOT/preflight2_${SLURM_JOB_ID}
mkdir -p "$OUT/gate" "$OUT/diag"
cp "$SUB"/case.json "$SUB"/*.py "$SUB"/job_preflight2.sh "$SUB"/gate/reference.json "$OUT/"
cp "$SUB"/gate/reference.json "$OUT/gate/"
sha256sum "$ROOT/build_v2/parties" "$ROOT/build_diag/parties" > "$OUT/binaries.sha256"

run() {  # run <dir> <ranks> <binary>
    ( cd "$1" && cp "$3" ./parties && \
      if srun --exact -n "$2" ./parties > run.log 2>&1; then echo "$(basename "$1") ok"; \
      else echo "$(basename "$1") FAILED"; fi >> "$OUT/status.txt" ) &
}
for v in G2_yang_n24 G2_yang_n32; do
    d=$OUT/gate/$v; mkdir -p "$d"; cp "$SUB/gate/$v/parties.inp" "$d/"
    printf '0\n' > "$d/p_mobile.inp"; printf '0\n' > "$d/p_fixed.inp"; printf '0\n' > "$d/stop.inp"
    run "$d" 8 "$ROOT/build_v2/parties"
done
for v in D1_base D2_gap3 D3_nograin D4_nobuoy D5_Re81 D6_yang D7_nomelt; do
    d=$OUT/diag/$v; mkdir -p "$d"; cp "$SUB/diag/$v/parties.inp" "$SUB/diag/$v/p_mobile.inp" "$d/"
    printf '0\n' > "$d/p_fixed.inp"; printf '0\n' > "$d/stop.inp"
    run "$d" 16 "$ROOT/build_diag/parties"
done
wait
echo "all steps finished" >> "$OUT/status.txt"
