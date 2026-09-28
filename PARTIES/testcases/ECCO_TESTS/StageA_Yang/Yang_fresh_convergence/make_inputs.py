#!/usr/bin/env python3
"""Generate the freshwater (Sm = 0) grid/interface-width convergence ladder.

Roadmap section A.5.2.  Every input is the completed Yang Sm=0 production input
(``Yang_production/parties.inp`` with cbd2 = cbd5 = 0) with ONLY the grid and
the four grid-tied numbers changed:

    dx = 1/N,  Cn = melt_band_eps = 0.75 dx,  Pe_CH = 0.9/Cn,
    zmax = twod_slab_thickness = dx,          dt scaled with dx.

Everything physical (Re, Pe, EOS, stefan, liquidus, darcy_tau, BCs, geometry,
ice slab position) is byte-identical across the ladder and identical to the
completed 1440^2 run, so the only variable is resolution.

Cn = 0.75 dx means the diffuse band narrows with the mesh: this is a
diffuse-interface convergence study, not just a mesh study.

GOTCHA: GRID_UNIFORM aborts unless dz matches dx to 1e-12 relative, so zmax is
written at full double precision (Python repr), never rounded.
"""
from pathlib import Path

HERE = Path(__file__).resolve().parent

# N: (ranks, wall-clock hours, nodes, partition)  -- Anvil, 128 cores/node.
# Rank counts factor cleanly into NPX x NPY x 1 (the code auto-factorizes):
#   360/64 = 8x8 -> 45x45 per rank; 512/128 = 16x8 -> 32x64; 720/256 = 16x16 -> 45x45.
# 'shared' is charged per core-hour, 'wholenode' per node-hour, so the two small
# grids go to 'shared' and only N720 takes whole nodes.
LADDER = {
    # Rungs BELOW the original ladder, added 2026-07-27.  Yang's base grid for
    # velocity and temperature in this case is 288^2 (their 5x refinement to
    # 1440^2 applies only to salinity and the phase field -- and the Sm=0 case
    # has NO salinity).  The ladder extrapolates to t_half ~ 102 at N=288
    # against Yang's 107.42, so these two rungs test whether the entire
    # "1.78x too fast" gap is PARTIES resolving heat delivery that Yang's base
    # grid does not.  288 = 2^5*3^2 and 216 = 2^3*3^3, both /8 exactly.
    216:  (64,   2, 1, "shared"),
    288:  (64,   3, 1, "shared"),
    360:  (64,   4, 1, "shared"),
    512:  (128,  8, 1, "shared"),
    # wholenode was backlogged ~2 days on 2026-07-27; 'shared' at 128 ranks
    # costs the same core-hours at half the parallelism and schedules now.
    # 720/16 = 45, 720/8 = 90, so 128 = 16x8 decomposes cleanly.
    720:  (128, 16, 1, "shared"),
}
TIME_MAX = 140.0        # past both t_half values (PARTIES 60.3, Yang 107.4)
OUT_DT = 2.0

TEMPLATE = """[geometry]
# A.5.2 freshwater convergence ladder, N = {N}.  Identical to the completed
# Sm=0 run (Yang A.4 geometry, s == 0) except for the grid and the four
# grid-tied numbers.  Yang 2-D reference for this case: t_half = 107.42,
# V(200)/V0 = 0.2483 (figure 3(a), dSv = 0 series).
xmin = 0.0
xmax = 1.0
ymin = 0.0
ymax = 1.0
zmin = 0.0
zmax = {dx!r}
twod_slab_thickness = {dx!r}
# = 1/{N}; keep FULL double precision or GRID_UNIFORM aborts at startup.

[grid]
NXM = {N}
NYM = {N}
NZM = 1
ImportGridFromFile = 0

[simulation]
time_max = {tmax}
output_time_interval = {odt}
output_time_interval_2d = {odt}
default_dt = {ddt!r}
max_dt = {mdt!r}
constant_dt = 0
cfl = 0.3
ghost_nodes = 3
resume = 0
freeze_velocity = 0

[flow]
Re = 1000.0
vel_init_type = 0
dp_dx = 0.0
ubulk_target = 0.0
startup_time = 0.0
startup_init = 0
startup_velocity = {{0.0, 0.0, 0.0}}

[vof]
We = 1.0e30
rho1 = 1.0
rho2 = 1.0
mu1 = 1.0
mu2 = 1.0
init_type = 27
vof_slab_x0 = 0.9
Cn = {cn!r}
# = 0.75*dx, the calibrated PARTIES mapping used by every passed gate
Pe_CH = {pech!r}
# = 0.9/Cn
darcy_tau = 1.0e-4
weno_order = 5
ch_iter_max = 3000
ch_tol = 1.0e-11
contact_angle_deg = 90.0

[phase_change]
stefan = 0.25
T_melt = 0.0
liquidus_slope = 0.014
melt_band_eps = {cn!r}

[eos]
eos_q = 2.0
eos_betaT = 1.0
eos_betaS = 1.75
eos_Tmd0 = 0.2
eos_Tmd_slope = -0.0625

[conc]
# NConc = 2 with s == 0 is kept (not NConc = 1) so this ladder is exactly
# comparable to the completed 1440^2 Sm=0 run, which carried the inert
# salt field.
NConc = 2
Nbin = 1000
Pe = {{1.0e4, 1.0e6}}
richardson = {{0.0, 0.0}}
V_s0 = {{0.0, 0.0}}
vol_heat_cap_s = {{1.0, 1.0}}
conductivity_s = {{1.0, 1.0}}
kappa_ice_ratio = {{1.0, 0.0}}
conc_init_type = {{30, 31}}
cbd0 = 1.0
cbd1 = 0.0
theta_ice = 0.0
cbd2 = 0.0
cbd5 = 0.0
# s_top = s_bot = 0: true freshwater, Sm = 0 and dSv = 0
BC_AW = {{1.0, 1.0}}
BC_BW = {{0.0, 0.0}}
BC_CW = {{0.0, 0.0}}
BC_AE = {{1.0, 1.0}}
BC_BE = {{0.0, 0.0}}
BC_CE = {{0.0, 0.0}}
BC_AN = {{1.0, 1.0}}
BC_BN = {{0.0, 0.0}}
BC_CN = {{0.0, 0.0}}
BC_AS = {{1.0, 1.0}}
BC_BS = {{0.0, 0.0}}
BC_CS = {{0.0, 0.0}}
BC_AF = {{1.0, 1.0}}
BC_BF = {{0.0, 0.0}}
BC_CF = {{0.0, 0.0}}
BC_AB = {{1.0, 1.0}}
BC_BB = {{0.0, 0.0}}
BC_CB = {{0.0, 0.0}}

[particle]
rho_s = 1.0
rho_prImp = 1.0
grav = {{0.0, -1.0, 0.0}}
N_forcing_loops = 1
N_heating_loops = 0
N_impheating_loops = 0

[lsolve]
P_CG_ETOL = 1.0e-10
P_CG_MAXIT = 2500
CG_ETOL = 1.0e-7
CG_MAXIT = 2500

[output]
ave_height_output = 0
front_location_output = 0
front_speed_output = 0
susp_mass_output = 0
sedim_rate_output = 0
energies_output = 0
shear_stress_output = 0
conc_output_integral = {{0, 0}}
conc_output_Q_int = {{0, 0}}
output_vfc = 0
output_vfu = 0
output_vfv = 0
output_vfw = 0
avg_dir = 1
nusselt_switch = 0
turbulent_terms_switch = 0
post_processing_switch = 0
near_wall_slice_switch = 0
slice_axis = 2
slice_half = 1
slice_pos = 1
slice_p = 0
center_two_particles = 0
"""

JOB = """#!/bin/bash -l
#SBATCH --job-name=YangFresh{N}
#SBATCH --output=fresh{N}.out
#SBATCH --error=fresh{N}.err
#SBATCH --account=phy250167
#SBATCH --partition={part}
#SBATCH --nodes={nodes}
#SBATCH --ntasks={ranks}
#SBATCH --time={hours}:00:00
#SBATCH --mail-user=mjalabert@ucsb.edu
#SBATCH --mail-type=END,FAIL

# Roadmap A.5.2 freshwater convergence ladder, N = {N}, t -> {tmax}.
# Compare t_half against the completed 1440^2 point (60.287) and the digitized
# Yang 2-D dSv=0 reference (107.42).

set -euo pipefail

# Anvil toolchain.  openmpi/4.1.6 is what the binary actually links against
# (check with `ldd parties`); HDF5 and HYPRE come from ~/Software/PARTIES_Libs
# via the binary's RPATH, so they need no module.  Do NOT use the gnu14/openmpi5
# names -- they do not exist on this machine and `set -e` will kill the job at
# the module load.
module load gcc/11.2.0 openmpi/4.1.6 fftw/3.3.8

PARTIES_BIN="${{PARTIES_BIN:-$HOME/PARTIES/PARTIES/parties}}"
if [ ! -x ./parties ]; then
    cp "$PARTIES_BIN" ./parties
fi

echo "Yang freshwater N={N} on $SLURM_NTASKS ranks, started $(date)"
srun ./parties > run.log 2>&1
echo "finished $(date)"
"""


def main():
    for N, (ranks, hours, nodes, part) in LADDER.items():
        d = HERE / f"N{N}"
        d.mkdir(exist_ok=True)
        dx = 1.0 / N
        cn = 0.75 / N
        scale = 1440.0 / N
        (d / "parties.inp").write_text(TEMPLATE.format(
            N=N, dx=dx, cn=cn, pech=0.9 / cn,
            tmax=TIME_MAX, odt=OUT_DT,
            ddt=1.0e-4 * scale, mdt=5.0e-4 * scale))
        jp = d / "jobscript.sh"
        jp.write_text(JOB.format(N=N, nodes=nodes, ranks=ranks, part=part,
                                 hours=hours, tmax=TIME_MAX))
        jp.chmod(0o755)
        # PARTIES reads stop.inp every iteration and aborts if it is absent
        # (Temporal_int.c:495).  Ship it with the deck.
        (d / "stop.inp").write_text("0\n")
        print(f"N={N:5d}  dx={dx:.10e}  Cn={cn:.10e}  Pe_CH={0.9 / cn:8.1f}  "
              f"dt0={1.0e-4 * scale:.3e}  ranks={ranks}  -> {d}")


if __name__ == "__main__":
    main()
