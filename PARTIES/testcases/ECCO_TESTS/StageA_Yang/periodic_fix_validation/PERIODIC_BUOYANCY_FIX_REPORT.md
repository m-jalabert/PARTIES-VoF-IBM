# Periodic buoyancy endpoint fix validation

## Change

`Velocity_add_buoyancy_2_RHS` now gives each face-velocity component the same
periodic endpoint extension already used by its normal momentum RHS.  The cell
centred EOS construction bounds are unchanged.  This removes the missing wrap
face that reduced `w` forcing to 3/4 on the four-cell periodic z grids.

## Validation

The two EOS cases were resumed from their original `t=12` states and relaxed
to `t=20` with the corrected executable.  The 4-dx penalization case was
resumed from its original `t=4` state and relaxed to `t=6.4` using `dt=0.002`.
Original HDF5 results were not modified.

| Case | Scale before | Scale after | Nominal relative error | Last-output change |
|---|---:|---:|---:|---:|
| Thermal EOS | 0.74999418 | 0.99986457 | 1.41e-4 | 6.05e-5 |
| Salty EOS | 0.74999768 | 0.99991996 | 8.39e-5 | 3.64e-5 |
| ICE penalization, 4 dx | 0.75000453 | 0.99997370 | 2.72e-5 | 3.22e-5 |

The thermal reconstructed EOS fit is `betaT=1.0002443`, `Tmd0=0.4000146`,
and `q=2.0004528`.  The salty reconstructed polynomial coefficients are
`c0=-0.3600050`, `c1=1.1499418`, and `c2=-0.2499509`, matching the nominal
`(-0.36, 1.15, -0.25)` coefficients.

The corrected 4-dx penalization profile matches the nominal diffuse Brinkman
BVP with scaled residual `1.36e-6`.  Its fitted decay length is
`0.0156655750`, versus discrete theory `0.0156655114`, a relative error of
`4.06e-6`.

## Conclusion

The periodic endpoint omission was the complete cause of the 0.75 amplitude in
both EOS and ICE-penalization gates.  With the correction:

- the nonlinear thermal/salinity EOS passes its absolute forcing test;
- the 4-dx ICE penalization case passes both absolute amplitude and decay-length
  tests;
- the previous conclusion that the Darcy coefficient itself was correct is
  confirmed.

The roadmap and earlier gate reports were not edited.
