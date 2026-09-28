#!/usr/bin/env python3
"""Salty (Sm=5, dSv=5) resolution test -- does the freshwater story transfer?

QUESTION
--------
The freshwater case closed because Yang's Sm=0 run has NO salinity, so their 5x
refinement bought nothing and 288^2 set their heat transport; PARTIES at 288^2
reproduced their t_half to 4.7%.

The salty case is NOT the same: Yang refines S and phi to 1440^2 there and only
u,T stay on 288^2.  PARTIES is single-grid, so neither of our options is a
faithful reproduction:

    PARTIES 288^2  -> matches their u,T resolution, but S and phi are 5x
                      COARSER than theirs (salt marginally resolved)
    PARTIES 1440^2 -> matches their S,phi, but u,T are 5x FINER (this is the
                      completed production run, t_half = 179.85)

So this is a diagnostic, not a reproduction.  Two outcomes, both informative:

  * lands near Yang's ~216  -> heat-transport resolution dominates the salty
    case too, the freshwater story transfers, and Stage A closes.
  * overshoots to ~300      -> the dual-grid SPLIT matters: coarsening S and phi
    along with u,T is not equivalent to Yang coarsening u,T alone.  Naive
    scaling of the freshwater ladder (t_half(288)/t_half(1440) = 1.698) predicts
    179.85 * 1.698 = 305, so this branch is the a-priori expectation.

Targets: Yang Sm=5,dSv=5 -> t_half ~ 216 (fig 3(b) f/f0 = 0.497 with
f0 = 1/107.42).  PARTIES 1440^2 -> 179.85.

RISK TO WATCH
-------------
Pe_S = 1e6 on a coarse mesh with 2nd-order CENTRAL scalar advection can produce
dispersive wiggles and out-of-range salinity.  The measured liquid-side salt
recovery length at 1440^2 was 27.7-98.4 cells, i.e. only ~5.5-19.7 cells at
288^2.  Check salt_min/salt_max and enthalpy drift before believing any
t_half from these runs -- the Step-C failure showed a melt rate can look clean
and decisive while the run is quietly non-conservative.
"""
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = HERE / "../Yang_production/parties.inp"

# N: (ranks, hours, nodes, partition).  288 = 2^5*3^2, 432 = 2^4*3^3; both /8 and /16.
LADDER = {
    288: (64, 6, 1, "shared"),
    432: (128, 10, 1, "shared"),
}
TIME_MAX = 350.0     # coarse rungs melt SLOWER; t_half may reach ~300
OUT_DT = 2.0

JOB = """#!/bin/bash -l
#SBATCH --job-name=YangSalt{N}
#SBATCH --output=salt{N}.out
#SBATCH --error=salt{N}.err
#SBATCH --account=phy250167
#SBATCH --partition={part}
#SBATCH --nodes={nodes}
#SBATCH --ntasks={ranks}
#SBATCH --time={hours}:00:00
#SBATCH --mail-user=mjalabert@ucsb.edu
#SBATCH --mail-type=END,FAIL

# Salty Sm=5 dSv=5 at N={N}, t -> {tmax}.  Targets: Yang ~216, PARTIES 1440^2 = 179.85.

set -euo pipefail
module load gcc/11.2.0 openmpi/4.1.6 fftw/3.3.8

PARTIES_BIN="${{PARTIES_BIN:-$HOME/PARTIES/PARTIES/parties}}"
if [ ! -x ./parties ]; then cp "$PARTIES_BIN" ./parties; fi

echo "Yang salty N={N} on $SLURM_NTASKS ranks, started $(date)"
srun ./parties > run.log 2>&1
echo "finished $(date)"
"""


def main():
    src = SRC.read_text().splitlines()
    for N, (ranks, hours, nodes, part) in LADDER.items():
        d = HERE / f"N{N}"
        d.mkdir(exist_ok=True)
        dx = 1.0 / N
        cn = 0.75 / N
        scale = 1440.0 / N
        out = []
        for line in src:
            s = line.strip()
            if s.startswith("zmax ="):
                out.append(f"zmax = {dx!r}")
            elif s.startswith("twod_slab_thickness ="):
                out.append(f"twod_slab_thickness = {dx!r}")
            elif s.startswith("NXM ="):
                out.append(f"NXM = {N}")
            elif s.startswith("NYM ="):
                out.append(f"NYM = {N}")
            elif s.startswith("Cn ="):
                out.append(f"Cn = {cn!r}")
            elif s.startswith("Pe_CH ="):
                out.append(f"Pe_CH = {0.9 / cn!r}")
            elif s.startswith("melt_band_eps ="):
                out.append(f"melt_band_eps = {cn!r}")
            elif s.startswith("time_max ="):
                out.append(f"time_max = {TIME_MAX}")
            elif s.startswith("output_time_interval ="):
                out.append(f"output_time_interval = {OUT_DT}")
            elif s.startswith("output_time_interval_2d ="):
                out.append(f"output_time_interval_2d = {OUT_DT}")
            elif s.startswith("default_dt ="):
                out.append(f"default_dt = {1.0e-4 * scale!r}")
            elif s.startswith("max_dt ="):
                out.append(f"max_dt = {5.0e-4 * scale!r}")
            else:
                out.append(line)
        hdr = (f"# Salty resolution test, N = {N}.  Identical to the completed 1440^2\n"
               f"# Sm=5 dSv=5 production input except the grid and the four grid-tied\n"
               f"# numbers (Cn = melt_band_eps = 0.75dx, Pe_CH = 0.9/Cn, zmax = dx,\n"
               f"# dt ~ dx) and time_max.  cbd2=0.5 / cbd5=1.5 keep the stratification.\n"
               f"# Cn = 0.75dx is a TUNED mapping -- do not vary it independently\n"
               f"# (see roadmap A.5.5 item 3, Step C failure).\n")
        (d / "parties.inp").write_text(hdr + "\n".join(out) + "\n")
        jp = d / "jobscript.sh"
        jp.write_text(JOB.format(N=N, nodes=nodes, ranks=ranks, part=part,
                                 hours=hours, tmax=TIME_MAX))
        jp.chmod(0o755)
        (d / "stop.inp").write_text("0\n")
        print(f"N={N:5d} dx={dx:.6e} Cn={cn:.6e} Pe_CH={0.9 / cn:7.1f} "
              f"dt0={1.0e-4 * scale:.3e} ranks={ranks} -> {d}")


if __name__ == "__main__":
    main()
