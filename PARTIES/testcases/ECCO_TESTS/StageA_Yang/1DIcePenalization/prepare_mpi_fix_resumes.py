#!/usr/bin/env python3
"""Set dt=0.002 in copied restart files for the 1-vs-4 rank check."""

from pathlib import Path
import h5py

ROOT = Path(__file__).resolve().parent
for run in ("tau_4dx_periodic_fix_mpi1", "tau_4dx_periodic_fix_mpi4"):
    filename = ROOT/run/"Resume.h5"
    with h5py.File(filename, "r+") as h5:
        h5["dt"][...] = 2.0e-3
        h5["dt_old"][...] = 2.0e-3
    print(filename)
