"""Tests for `equations.py` helper functions and stress calculations."""
import os
import unittest

import numpy as np

from stress_analysis import data_ingestion as data_in
from stress_analysis import equations as eqns


class TestInterpolateX(unittest.TestCase):
    """Tests for the `interpolate_xu_to_xc` method."""

    def test_simple_field(self):
        """Test the interpolation of a simple field along the `x`-axis."""
        test_field = np.array(
            [
                [
                    [0, 1],
                    [0, 1],
                ],
                [
                    [0, 1],
                    [0, 1],
                ],
            ]
        )
        result = eqns.interpolate_xu_to_xc(test_field)
        expected = np.array(
            [
                [
                    [0.5],
                    [0.5],
                ],
                [
                    [0.5],
                    [0.5],
                ],
            ]
        )
        np.testing.assert_almost_equal(result, expected)


class TestInterpolateY(unittest.TestCase):
    """Tests for the `interpolate_yv_to_yc` method."""

    def test_simple_field(self):
        """Test the interpolation of a simple field along the `y`-axis."""
        test_field = np.array(
            [
                [
                    [1, 1],
                    [0, 0],
                ],
                [
                    [1, 1],
                    [0, 0],
                ],
            ]
        )
        result = eqns.interpolate_yv_to_yc(test_field)
        expected = np.array(
            [
                [
                    [0.5, 0.5],
                ],
                [
                    [0.5, 0.5],
                ],
            ]
        )
        np.testing.assert_almost_equal(result, expected)


class TestAddGhostCoords(unittest.TestCase):
    """Test `add_ghost_coords` method."""

    def test_simple_coords(self):
        """Test the addition of ghost nodes to a 1D coordinate array."""
        test_array = np.array([0, 1, 2, 3])
        h = 0.5
        result = eqns.add_ghost_coords(test_array, h)
        expected = np.array([-0.5, 0, 1, 2, 3, 3.5])
        np.testing.assert_almost_equal(result, expected)


class TestAddGhostNodes(unittest.TestCase):
    """Test `add_ghost_nodes` method."""

    def test_default_args(self):
        """Test without specifying bulk velocities."""
        test_field = np.array(
            [
                [
                    [0, 1],
                    [0, 1],
                ],
                [
                    [0, 1],
                    [0, 1],
                ],
            ]
        )
        result = eqns.add_ghost_nodes(test_field)
        expected = np.array(
            [
                [
                    [0, -1],
                    [0, 1],
                    [0, 1],
                    [0, -1],
                ],
                [
                    [0, -1],
                    [0, 1],
                    [0, 1],
                    [0, -1],
                ],
            ]
        )
        np.testing.assert_almost_equal(result, expected)

    def test_specified_u_bulk(self):
        """Test with specified bulk velocities."""
        test_field = np.array(
            [
                [
                    [0, 1],
                    [0, 1],
                ],
                [
                    [0, 1],
                    [0, 1],
                ],
            ]
        )
        u_bulk = 1
        result = eqns.add_ghost_nodes(test_field, top=u_bulk, bot=u_bulk)
        expected = np.array(
            [
                [
                    [2, 1],
                    [0, 1],
                    [0, 1],
                    [2, 1],
                ],
                [
                    [2, 1],
                    [0, 1],
                    [0, 1],
                    [2, 1],
                ],
            ]
        )
        np.testing.assert_almost_equal(result, expected)


class TestStressX(unittest.TestCase):
    """Integration test for the `stress_x` function."""

    def test_balance(self):
        """Test the equality of the stress balance in the x direction."""
        test_dir = os.path.dirname(os.path.realpath(__file__))
        data_dir = os.path.join(test_dir, "test_data", "rolling_particle")
        sim = data_in.Simulation(data_dir)
        stress_x = eqns.stress_x(sim, sim.fluid.indices[0])
        np.testing.assert_allclose(
            stress_x["total"], stress_x["external"], rtol=0.1
        )


class TestStressY(unittest.TestCase):
    """Integration test for the `stress_y` function."""

    def test_balance(self):
        """Test the equality of the stress balance in the y direction."""
        test_dir = os.path.dirname(os.path.realpath(__file__))
        data_dir = os.path.join(test_dir, "test_data", "rolling_particle")
        sim = data_in.Simulation(data_dir)
        stress_y = eqns.stress_y(sim, sim.fluid.indices[0])
        np.testing.assert_allclose(
            stress_y["total"], stress_y["external"], rtol=0.1
        )


class TestParticlePressure(unittest.TestCase):
    """Integration test for the `particle_pressure` function."""

    def test_particle_weight(self):
        """Test if the particle weight is correct for the example simulation."""
        test_dir = os.path.dirname(os.path.realpath(__file__))
        data_dir = os.path.join(test_dir, "test_data", "particle_on_bed")
        sim = data_in.Simulation(data_dir)
        p_pressure = eqns.particle_pressure(sim, sim.p_fixed.indices[0])
        grav = sim.grav[1]
        rho_s = sim.rho_s
        area = sim.Lx * sim.Lz
        expected_pressure = (
            (4 / 3 * np.pi * sim.p_fixed.R(sim.p_fixed.indices[0])[0, 0] ** 3)
            * (rho_s - 1)
            * grav
            / area
        )
        np.testing.assert_allclose(expected_pressure, p_pressure["bot"])
        np.testing.assert_allclose(0, p_pressure["top"])
        np.testing.assert_allclose(expected_pressure / 2, p_pressure["avg"])


if __name__ == "__main__":
    unittest.main()
