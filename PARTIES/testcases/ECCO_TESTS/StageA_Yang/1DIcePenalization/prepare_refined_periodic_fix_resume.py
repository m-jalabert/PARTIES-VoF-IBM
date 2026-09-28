#!/usr/bin/env python3
"""Set dt=0.002 in the copied t=4 penalization restart."""

from pathlib import Path
import h5py

filename = Path(__file__).resolve().parent / "tau_4dx_periodic_fix_refined/Resume.h5"
with h5py.File(filename, "r+") as h5:
    h5["dt"][...] = 2.0e-3
    h5["dt_old"][...] = 2.0e-3
print(filename)
