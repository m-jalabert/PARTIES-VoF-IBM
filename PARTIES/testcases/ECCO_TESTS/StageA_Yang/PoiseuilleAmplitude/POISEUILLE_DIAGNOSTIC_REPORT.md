# Fixed-pressure-gradient Poiseuille diagnostic

## Purpose

Determine whether the common measured amplitude of approximately 0.75 in the
EOS and ICE-penalization gates is caused by an over-strong momentum operator or
by the buoyancy-source path.

The diagnostic uses an all-liquid, equal-property diffuse-VOF field, with
concentration, Boussinesq forcing, nonlinear EOS, phase change, and ice
penalization disabled.  The flow is driven only by `dp_dx = -1` in a periodic
streamwise direction, with `Re = 10` and no-slip walls at `y = 0,1`.

## Results

- The 256-cell exact-profile-start run gives measured/discrete-exact amplitude
  `0.9999973425` at `t = 2`, with relative L-infinity profile error
  `2.75e-6` and reconstructed forcing `0.9999978427`.
- The independent 64-cell zero-start run reaches amplitude `0.9996276766` at
  `t = 8`.  A late-time exponential fit gives asymptotic amplitude
  `0.9999999999`.
- The fitted decay rate is `0.9867867627`, within `1.76e-4` relative error of
  the theoretical slowest viscous-mode rate `pi^2/Re = 0.9869604401`.
- Both transverse velocity components are exactly zero to stored precision,
  and streamwise/spanwise nonuniformity is at roundoff level.
- The phase indicator remains exactly one in both clean runs.

## Conclusion

The momentum diffusion/operator normalization is correct.  It does not produce
the factor 0.75.  The fixed pressure-gradient source balances the discrete
viscous operator at amplitude one.

Source inspection explains the prior factor.  The EOS and penalization gates
drive `w` using gravity in periodic `z`, with `NZM = 4`.  The normal `w` RHS
extends `k_end` for `ZPERIODIC`, but `Velocity_add_buoyancy_2_RHS` computes
`k_end = min(NZ-1, Ke)` without the periodic extension before looping over the
`w` faces.  Consequently one of four periodic velocity planes does not receive
the direct buoyancy source.  The projected mean forcing is therefore 3/4 of
nominal, matching the measured `0.7500` in both prior gates.

This means:

- the nonlinear EOS algebraic form is correct, but the gate's realized force
  amplitude is reduced by the periodic-loop-bound defect;
- the ICE penalization coefficient and decay-length scaling remain correct;
  its absolute velocity amplitude is low because its driving buoyancy follows
  the same defective loop, not because the Darcy term is too strong.

No solver-source correction was made as part of this diagnostic.
