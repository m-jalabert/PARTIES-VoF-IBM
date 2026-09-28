#!/usr/bin/env python3
"""Build a nonzero-time multicomponent Stefan seed from a PARTIES output.

The AFiD-MuRPhFi validation starts from the continuous similarity profiles at
nonzero time.  This generator keeps the PARTIES HDF5 schema and replaces only
the phase, temperature, salinity, zero-flow state, and restart counters.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
import shutil

import h5py
import numpy as np

from analyze_gate import (
    DELTA,
    LE,
    PE_S,
    PE_T,
    THETA_ICE,
    THETA_INF,
    X0,
    reference_state,
)


CN = 0.00146484375
SEED_TIME = 0.25
DT = 1.0e-4
SEED_OUTPUT = 5


def similarity_fields(x: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    lam, s_i, theta_i = reference_state()
    front = X0 + 2.0 * lam * math.sqrt(SEED_TIME / PE_T)

    f = 0.5 * (1.0 - np.tanh((x - front) / (2.0 * math.sqrt(2.0) * CN)))

    erf = np.vectorize(math.erf)
    erfc = np.vectorize(math.erfc)
    eta_t = (x - X0) / (2.0 * math.sqrt(SEED_TIME / PE_T))
    eta_s = (x - X0) / (2.0 * math.sqrt(SEED_TIME / PE_S))

    theta_liquid = THETA_INF + (theta_i - THETA_INF) * (1.0 + erf(eta_t)) / (
        1.0 + math.erf(lam)
    )
    theta_solid = THETA_ICE + (theta_i - THETA_ICE) * erfc(eta_t) / math.erfc(lam)
    theta = np.where(x <= front, theta_liquid, theta_solid)

    lam_s = lam * math.sqrt(LE)
    salt = 1.0 + (s_i - 1.0) * (1.0 + erf(eta_s)) / (
        1.0 + math.erf(lam_s)
    )

    return f, theta, salt


def broadcast(profile: np.ndarray, shape: tuple[int, ...]) -> np.ndarray:
    return np.broadcast_to(profile.reshape((1,) * (len(shape) - 1) + (-1,)), shape)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    source_data = args.source / "Data_0.h5"
    source_resume = args.source / "Resume_0.h5"
    args.output.mkdir(parents=True, exist_ok=False)
    data_path = args.output / f"Data_{SEED_OUTPUT}.h5"
    resume_path = args.output / "Resume.h5"
    shutil.copy2(source_data, data_path)
    shutil.copy2(source_resume, resume_path)

    with h5py.File(data_path, "r+") as h5:
        x = np.asarray(h5["grid/xc"], dtype=np.float64)
        f, theta, salt = similarity_fields(x)

        for name in ("VOF/F", "VOF/C_L"):
            h5[name][...] = broadcast(f, h5[name].shape)
        h5["VOF/C_G"][...] = broadcast(1.0 - f, h5["VOF/C_G"].shape)
        h5["VOF/C_S"][...] = 0.0
        h5["VOF/ch_rhs_n"][...] = 0.0
        h5["VOF/ch_rhs_nm1"][...] = 0.0
        h5["Conc/0"][...] = broadcast(theta, h5["Conc/0"].shape)
        h5["Conc/1"][...] = broadcast(salt, h5["Conc/1"].shape)
        for name in ("u", "v", "w", "p"):
            h5[name][...] = 0.0
        h5["time"][...] = SEED_TIME

    with h5py.File(resume_path, "r+") as h5:
        fshape = h5["Resume/C_L"].shape
        f, _, _ = similarity_fields(x)
        for name in ("Resume/F", "Resume/C_L"):
            h5[name][...] = broadcast(f, fshape)
        h5["Resume/C_S"][...] = 0.0
        h5["Resume/ch_rhs_n"][...] = 0.0
        h5["Resume/ch_rhs_nm1"][...] = 0.0
        h5["time"][...] = SEED_TIME
        h5["dt"][...] = DT
        h5["dt_old"][...] = DT
        h5["ntime"][...] = int(round(SEED_TIME / DT)) + 1
        h5["noutput"][...] = SEED_OUTPUT
        h5["output_time"][...] = SEED_TIME

    lam, s_i, theta_i = reference_state()
    print(
        f"seed={data_path} time={SEED_TIME} lambda={lam:.12g} "
        f"s_i={s_i:.12g} theta_i={theta_i:.12g} delta={DELTA:g}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
