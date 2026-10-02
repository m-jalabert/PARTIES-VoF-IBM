#!/bin/bash -l
# dev3 (debug queue): Yang-operator repair at resolved sediment (build_v3).
#   Y_orig  small box, Yang, build_v2 (expected to diverge, as P1b job 20973340)
#   Y_fix   same deck, build_v3 (repair: capacity C_L+C_S under VOF_IBM)
#   G2r     1-D Yang gate G2_yang_n24 with build_v3 (no grain: must match round 2)
#   NI_v2/NI_v3  default operator, build_v2 vs build_v3 (must be bit-identical)
#SBATCH --job-name=P1b_dev3
#SBATCH --account=phy250167
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --time=01:00:00
#SBATCH --output=/anvil/scratch/x-mjalabert/ECCO_StageB3/P1_Attempts/P1b_shelf/logs/dev3_%j.out
#SBATCH --error=/anvil/scratch/x-mjalabert/ECCO_StageB3/P1_Attempts/P1b_shelf/logs/dev3_%j.err
set -uo pipefail
module purge
module load gcc/11.2.0 openmpi/4.1.6 hdf5/1.10.7 fftw/3.3.8
ROOT=/anvil/scratch/x-mjalabert/ECCO_StageB3/P1_Attempts/P1b_shelf
SUB=$SLURM_SUBMIT_DIR
OUT=$ROOT/dev3_${SLURM_JOB_ID}
mkdir -p "$OUT"
cp "$SUB"/case.json "$SUB"/*.py "$SUB"/job_dev3.sh "$OUT/"
sha256sum "$ROOT"/build_v2/parties "$ROOT"/build_v3/parties > "$OUT/binaries.sha256"
run() {  # run <name> <deck dir> <ranks> <build>
    d=$OUT/$1; mkdir -p "$d"; cp "$2"/parties.inp "$d/"
    if [ -f "$2/p_mobile.inp" ]; then cp "$2/p_mobile.inp" "$d/"; else printf '0\n' > "$d/p_mobile.inp"; fi
    printf '0\n' > "$d/p_fixed.inp"; printf '0\n' > "$d/stop.inp"; cp "$ROOT/$4/parties" "$d/parties"
    ( cd "$d" && if srun --exact -n "$3" ./parties > run.log 2>&1; then echo "$1 ok"; else echo "$1 FAILED"; fi >> "$OUT/status.txt" ) &
}
run Y_orig "$SUB/dev3/yang" 48 build_v2
run Y_fix  "$SUB/dev3/yang" 48 build_v3
run G2r    "$SUB/gate/G2_yang_n24" 8 build_v3
run NI_v2  "$SUB/dev3/noninv" 12 build_v2
run NI_v3  "$SUB/dev3/noninv" 12 build_v3
wait
echo "all steps finished" >> "$OUT/status.txt"
