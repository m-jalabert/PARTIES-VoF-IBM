#!/usr/bin/env python3
"""Set dt=0.002 in copied t=4 restart files for corrected 2/8-dx runs."""

from pathlib import Path
import h5py

ROOT = Path(__file__).resolve().parent
for run in ("tau_2dx_periodic_fix_final", "tau_8dx_periodic_fix_final"):
    filename = ROOT / run / "Resume.h5"
    with h5py.File(filename, "r+") as h5:
        h5["dt"][...] = 2.0e-3
        h5["dt_old"][...] = 2.0e-3
    print(filename)
