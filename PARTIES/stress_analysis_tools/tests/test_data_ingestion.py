"""Tests for `data_ingestion.py` helper functions."""
import os
import tempfile
import unittest
import warnings

import numpy as np

from stress_analysis import data_ingestion as data_in


class TestFindParameter(unittest.TestCase):
    """Test the `find_parameter` method."""

    def test_find_parameter(self):
        """Test if numerical parameters are correctly read from a file."""
        # Write temporary parties.inp file
        param_name = "sample_param"
        param_value = 1.5
        sample_text = (
            f"# This is a comment\n{param_name} = {param_value}\nNYM = 1\n"
        )
        with tempfile.NamedTemporaryFile() as tmp:
            with open(tmp.name, "w", encoding="utf8") as file:
                file.write(sample_text)
            with open(tmp.name, "r", encoding="utf8") as file:
                with self.subTest():
                    result = data_in.find_parameter(file, param_name)
                    self.assertAlmostEqual(
                        param_value,
                        result,
                        6,
                        "Failure to read defined `parties.inp` parameters.",
                    )
                with self.subTest():
                    result = data_in.find_parameter(file, "non_param")
                    self.assertAlmostEqual(
                        0.0,
                        result,
                        6,
                        "Failed managing undefined `parties.inp` parameter.",
                    )


class TestFindNumber(unittest.TestCase):
    """Test the `find_number` method."""

    def setUp(self):
        warnings.simplefilter("ignore", category=UserWarning)

    def test_float(self):
        """Test if simple float is found."""
        line = "xmax = 1.0"
        expected = 1.0
        result = data_in.find_number(line)
        np.testing.assert_allclose(expected, result)

    def test_scientific_notation(self):
        """Test if number written in scientific notation is found."""
        line = "xmax = 1e-3"
        expected = 1e-3
        result = data_in.find_number(line)
        np.testing.assert_almost_equal(expected, result)

    def test_multiple_numbers(self):
        """Test if a vector of numbers is found and read correctly."""
        line = "grav = {12, 0.1, 14}"
        expected = np.array([12, 0.1, 14])
        result = data_in.find_number(line)
        np.testing.assert_allclose(result, expected)

    def test_missing(self):
        """Test if there is no number in the input string."""
        line = "xmax = "
        expected = 0.0
        result = data_in.find_number(line)
        np.testing.assert_almost_equal(expected, result)


class TestGetIndices(unittest.TestCase):
    """Test the `get_indices` method."""

    def test_for_existing_files(self):
        """Test for a directory filled with multiple files."""
        wildcard = "Data_*.h5"
        test_dir = os.path.dirname(os.path.realpath(__file__))
        data_dir = os.path.join(test_dir, "test_data", "rolling_particle")
        expected = np.array([5, 6, 7])
        result = data_in.get_indices(data_dir, wildcard)
        np.testing.assert_almost_equal(result, expected)

    def test_missing_files(self):
        """Test for behavior when no files match the wildcard."""
        wildcard = "NoFile_*.h5"
        test_dir = os.path.dirname(os.path.realpath(__file__))
        data_dir = os.path.join(test_dir, "test_data", "rolling_particle")
        expected = np.array([])
        result = data_in.get_indices(data_dir, wildcard)
        np.testing.assert_almost_equal(result, expected)


if __name__ == "__main__":
    unittest.main()
