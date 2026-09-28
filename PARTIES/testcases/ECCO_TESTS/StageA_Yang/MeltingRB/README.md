# Melting Rayleigh-Benard benchmark (Favier, Purseed & Duchemin, JFM 858, 2019)

Purpose: salt-free, convection-coupled melting benchmark against an
independent published reference, to discriminate the Yang A.4 melt-rate
deficit (see `docs/md_files/YANG_PRODUCTION_ANALYSIS_REPORT.md`): the
interface machinery is proven exact (gates + in-situ audits) and
tau/Pe_CH-insensitive, so the remaining suspects are the convective heat
delivery or the Yang reference itself. This case tests exactly the
convective melt coupling, with no salt.

Paper case D: `Ra = 1e7`, `Pr = 1`, `thetaM = 0.05`, `h0 = 0.05`,
aspect `lambda = 6`, horizontal periodic, no-slip walls, solid initially
near the melting temperature (paper: docs/papers/Favier2019_MeltingRB.pdf).

## PARTIES realization

- `Boundary.h`: `PERIODIC_NOSLIP_BOX` (x periodic, y no-slip walls, z
  periodic storage slab) + `TWOD_CARTESIAN` + the full Yang flag set
  unchanged. **Restore `PERIODIC_Z_NOSLIP_BOX` for Yang A.4 runs.**
- New IC types (this campaign): VOF `init_type = 28` (horizontal solid
  layer above liquid, `vof_slab_x0` = h0) and Conc `init_type = 33`
  (Favier eq. 4.2 piecewise-linear temperature, `T_melt` = thetaM,
  `cbd1` = perturbation amplitude in the liquid only; deterministic modes).
- Linear buoyancy via the (validated) nonlinear-EOS path degenerated to
  q = 1, betaT = 1, betaS = 0, Tmd0 = 0, Tmd_slope = 0  =>  b = -theta.
- Free-fall units: `Re = sqrt(Ra/Pr) = 3162.2777`, `Pe_T = 3162.2777`.
  Favier's diffusive-time results convert as t_diff = t_ff/Pe_T.
- `stefan` (= cp dT/L) = 1/St_Favier: run St1 (stefan=1), St0p1 (stefan=10).
- Grid 3072 x 512 (dx = 1/512, cubic), Cn = 0.75 dx, Pe_CH = 0.9/Cn,
  melt_band_eps = Cn, darcy_tau = 1e-4, mu2 = 1 (all production-consistent).

## Quantitative pass criteria (paper Fig. 7, same Ra)

**CORRECTED 2026-07-17** (the original criterion 2 divided by St instead of
(1/2+St); paper eq. (48)-(51) include the sensible-heat storage of the newly
molten fluid, mean temperature 1/2):

1. Diffusive phase: h(t) follows the exact 1-D two-phase Stefan solution
   (integrate it numerically for the eq. 4.2 IC; asymptotically
   h^2 ~ h0^2 + 4 lam^2 t_diff with lam(St=1, thetaM=0.05) = 0.608)
   until Ra_e = Ra (1-thetaM) h^3 = 1708  (h_c = 0.0564).
2. Post-onset: h grows linearly; melting velocity (diffusive units)
   (1/2 + St) hdot = gamma Ra^{1/3} (1-thetaM)^{4/3} with gamma = 0.115
   => hdot(St=1) = 23.14/1.5 = 15.4, hdot(St=0.1) = 23.14/0.6 = 38.6:
   the St-dependence of Fig. 7(b) (ratio 2.5, saturating at low St, NOT 10x).
3. Nu(Ra_e) consistent with gamma Ra_e^{1/3} (their Fig. 10).

## Result (St=1 run, 2026-07-17, Anvil job 19301354)

Run `/anvil/scratch/x-mjalabert/Melting_RB_07162026`, t -> 120, clean.
Diffusive phase: 1.7% (PASS). Onset: consistent with h_c (PASS). Interface
energy conversion (q_wall vs (1/2+St) Pe hdot): closes to 0.4-1.7% (PASS -
kernel exact). Melt plateau: 18.3-19.8 vs 23.14 (ratio 0.79-0.85, MISS);
Nu amplitude gamma_eff = 0.0955 +/- 0.009 vs 0.115 (-17%), exponent 0.348
(shape PASS). Morphology: scallops -> 10 topography-locked cells (PASS).
Full analysis: `docs/md_files/NASA-JPL/MELTING_RB_ANALYSIS_REPORT.md`.

Gate script: `analyze_melting_rb.py` (writes per-run timeseries CSVs and
prints diffusive-phase error, onset height, post-onset hdot vs target).

Runs live in `/bigscratch/mjalabert314/MeltingRB/{St1,St0p1}`.

## Verdict logic

- PARTIES matches Fig. 7 => convective thermal melt coupling correct =>
  Yang deficit localizes to the salt-side flow structure or the reference;
  next: fresh (Sm=0) Yang run for the paper's f/f0 ratio test.
- PARTIES ~30% slow here too => minimal salt-free reproduction of the
  Yang deficit; debug convective melt coupling at matched parameters
  (Allen-Cahn vs Cahn-Hilliard comparison becomes concrete).
