#!/bin/bash -l
#SBATCH --job-name=SAgateB14
#SBATCH --output=gate_b14.out
#SBATCH --error=gate_b14.err
#SBATCH --account=phy250167
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=00:30:00
#
# Roadmap B.1.0 mandated re-run: "Re-run this gate whenever a new ECCO change
# lands on a Stage-A-reachable file (Conc.c, VOF_DIFFUSE.c, Initial_Conditions.c,
# Velocity.c)."  B.1.4 touches VOF_DIFFUSE.c, so this is required, not optional.
#
# Binary parties.stageA_b14: current source, ECCO flags OFF.  Must reproduce the
# recorded Le=100 salt gate to every printed digit:
#     lambda = 0.21844955   rel.err 6.039e-03   profile L_inf 1.441e-02
#     salt drift 2.764e-03  max speed 0
#
# Stage-A build has LAG_PARTICLE_RESOLVED off, so NO p_*.inp are needed.
set -euo pipefail
set +u; module purge; module load gcc/11.2.0 openmpi/4.1.6 hdf5/1.10.7 fftw/3.3.8; set -u
C="${SLURM_SUBMIT_DIR}"
# Run output goes to SCRATCH, never $HOME (25 GB quota).  See roadmap B.0.1.
d="/anvil/scratch/x-mjalabert/ECCO_StageB/StageA_noimpact_gate_b14_${SLURM_JOB_ID}"
mkdir -p "$d"
cp "$C/parties.inp" "$C/stop.inp" "$d/"
cp /home/x-mjalabert/PARTIES/PARTIES/parties.stageA_b14 "$d/parties"
cd "$d" && srun ./parties > run.log 2>&1
echo "stageA no-impact gate (B.1.4 recheck) done"
