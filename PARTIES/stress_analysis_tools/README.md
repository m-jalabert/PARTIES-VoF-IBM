# Stress Analysis

Perform stress analysis on PARTIES simulation data.

The primary purpose of the module is to provide a user friendly way to ingest
PARTIES data and calculate x and y stress profiles by implementing equations
7.12 and 7.14 in Edward Biegert's PhD disertation (2018). The module can also
be used to calculate other quantities, such as particle pressure.

Dependencies
------------

The simple way to run the module is to create a conda environment using the
provided `environment.yml` file. To do so use:
`conda env create -f environment.yml`

Optionally one can create their own environment and install the dependencies
in the `environment.yml` file manually.

Requires python 3.8+

Usage
-----

See `examples` for a jupyter notebook demonstrating typical usage.

Testing
-------

Tests can be run with:
`python -m unittest tests/test_data_ingestion.py`
`python -m unittest tests/test_equation.py`