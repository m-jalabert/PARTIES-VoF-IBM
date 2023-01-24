"""Data structures and methods to handle parties h5 data."""
import glob
import os
import re
from typing import TextIO

import h5py
import numpy as np


class Simulation:
    """PARTIES simulation parameters, geometry, and data.

    PARTIES simulation parameters and geometry are read from the `parties.inp`
    file found alongside the data output. Fluid and particle data are stored in
    respective sub-structures.
    """

    def __init__(self, datapath: str) -> None:
        self.datapath = datapath

        self.Lx = self.get_input_parameter("xmax")
        self.Ly = self.get_input_parameter("ymax")
        self.Lz = self.get_input_parameter("zmax")
        self.Nx = self.get_input_parameter("NXM").astype("int")
        self.Ny = self.get_input_parameter("NYM").astype("int")
        self.Nz = self.get_input_parameter("NZM").astype("int")
        self.Re = self.get_input_parameter("Re")
        self.ubulk = self.get_input_parameter("ubulk_target")
        self.dp_dx = -self.get_input_parameter("dp_dx")
        self.rho_s = self.get_input_parameter("rho_s")
        self.grav = self.get_input_parameter("grav")
        self.h = self.Ly / self.Ny

        self.fluid = FluidData(datapath)
        self.p_fixed = ParticleData(datapath, "fixed")
        self.p_mobile = ParticleData(datapath, "mobile")

    def get_input_parameter(self, quantity: str) -> float:
        """Find input parameters in a `parties.inp` file.

        Args:
            quantity: The desired input parameter as written in `parties.inp`.

        Returns:
            The value of the quantity if it exists. `0.0` if not.
        """
        parties_inp_filepath = os.path.join(self.datapath, "parties.inp")
        with open(os.path.join(parties_inp_filepath), encoding="utf8") as file:
            return find_parameter(file, quantity)


class FluidData:
    """Container and methods for handilng Data_*.h5 parties output."""

    def __init__(self, datapath: str) -> None:
        self.datapath = datapath
        self.indices = get_indices(datapath, "Data_*.h5")
        if self.indices.size > 0:
            filename = os.path.join(datapath, f"Data_{self.indices[0]}.h5")
            with h5py.File(filename) as file:
                self.xc = file["grid"]["xc"][:-1]
                self.yc = file["grid"]["yc"][:-1]
                self.zc = file["grid"]["zc"][:-1]
                self.xu = file["grid"]["xu"][:]
                self.yv = file["grid"]["yv"][:]
                self.zw = file["grid"]["zw"][:]

    def u(self, idx: int) -> np.ndarray:
        """Read x fluid velocity field definded at cell edges in x direction."""
        return self.read_field(idx, "u", "u")

    def v(self, idx: int) -> np.ndarray:
        """Read y fluid velocity field definded at cell edges in y direction."""
        return self.read_field(idx, "v", "v")

    def w(self, idx: int) -> np.ndarray:
        """Read z fluid velocity field definded at cell edges in z direction."""
        return self.read_field(idx, "w", "w")

    def p(self, idx: int) -> np.ndarray:
        """Read pressure field definded at cell centers."""
        return self.read_field(idx, "p", "c")

    def fx_ibm(self, idx: int) -> np.ndarray:
        """Read local ibm forcing in the x-direction, defined at x cell edges."""
        return self.read_field(idx, "fx_IBM", "u")

    def fy_ibm(self, idx: int) -> np.ndarray:
        """Read local ibm forcing in the y-direction, defined at y cell edges."""
        return self.read_field(idx, "fy_IBM", "v")

    def vfu(self, idx: int) -> np.ndarray:
        """Read particle indicator function defined at x cell edges."""
        return self.read_field(idx, "vfu", "u")

    def vfv(self, idx: int) -> np.ndarray:
        """Read particle indicator function defined at y cell edges."""
        return self.read_field(idx, "vfv", "v")

    def vfw(self, idx: int) -> np.ndarray:
        """Read particle indicator function defined at z cell edges."""
        return self.read_field(idx, "vfw", "w")

    def vfc(self, idx: int) -> np.ndarray:
        """Read cell centered particle indicator function."""
        return self.read_field(idx, "vfc", "c")

    def read_field(self, idx: int, field: str, grid: str) -> np.ndarray:
        """Read in desired fluid data field and remove ghost nodes.

        Args:
            idx: A h5 file output index.
            field: Whcih h5 attribute to read from file.
            grid: Which grid the field is defined.
        Returns:
            The desired field without ghost nodes.
        """
        filename = f"Data_{idx}.h5"
        with h5py.File(os.path.join(self.datapath, filename)) as file:
            if grid == "u":
                return file[field][:-1, :-1, :]
            if grid == "v":
                return file[field][:-1, :, :-1]
            if grid == "w":
                return file[field][:, :-1, :-1]
            if grid == "c":
                return file[field][:-1, :-1, :-1]
            return file[field][:]


class ParticleData:
    """Fixed or mobile particle data and methods.

    Note:
        Particle objects are made even if there are no `fixed` or `mobile`
        particles in the simulation. This can be worked on in a future update.

    Attributes:
        datapath: A string specifying the directory where `Particle_*.h5` files
            are stored.
        particle_type: `fixed` or `mobile`, a string specifying particle type.
    """

    def __init__(self, datapath, particle_type):
        self.datapath = datapath
        self.indices = get_indices(datapath, "Particle_*.h5")
        self.particle_type = particle_type
        if self.indices.size > 0:
            try:
                self.count = len(self.X(self.indices[0]))
            except KeyError:
                self.count = 0
                self.indices = np.array([])

    def X(self, idx: int) -> np.ndarray:
        """Read particle position vectors."""
        return self.read_data(idx, "X")

    def U(self, idx: int) -> np.ndarray:
        """Read particle velocity vectors."""
        return self.read_data(idx, "U")

    def R(self, idx: int) -> np.ndarray:
        """Read particle radii."""
        return self.read_data(idx, "R")

    def F_coll(self, idx: int) -> np.ndarray:
        """Read particle collision force vectors."""
        return self.read_data(idx, "F_coll")

    def read_data(self, idx: int, field: str) -> np.ndarray:
        """Read in desired particle data field.

        Args:
            idx: A  h5 file output index.
            field: Whcih h5 attribute to read from file.
        Returns:
            The desired particle data.
        """
        filename = f"Particle_{idx}.h5"
        with h5py.File(os.path.join(self.datapath, filename)) as file:
            return file[self.particle_type][field][:]


def get_indices(datapath: str, wildcard: str) -> np.ndarray:
    """Get h5 file indices for specified wildcard.

    Example:
        get_indices('Data_*.h5')
        get_indices('Particle_*.h5')

    Args:
        datapath: The directory containing parties h5 output data.
        wildcard: The desired h5 file wildcards. See examples.

    Returns:
        The corresponding h5 file indices found in the given directory.
    """
    filepaths = glob.glob(os.path.join(datapath, wildcard))
    filenames = [os.path.basename(filepath)[:-3] for filepath in filepaths]
    indices = np.array([int(find_number(filename)) for filename in filenames])
    indices.sort()
    return indices


def find_parameter(file: TextIO, quantity: str) -> np.ndarray:
    """In a file, find a number on the same line as a desired word or phrase."""
    lines = [ln for ln in file.read().split("\n") if quantity in ln]
    if lines:
        return find_number(lines[0])
    return np.zeros(1)


def find_number(line: str) -> np.ndarray:
    """Extract the first number in a string."""
    regex = r"-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?"
    return np.asarray(re.findall(regex, line), dtype=float)
