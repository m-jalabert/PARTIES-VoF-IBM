# Documentation for Mordant-Testcase

## Files of the testcase
Boundary_Mordant.h (corresponds to the Boundary.h-File), p_fixed.inp, p_mobile.inp, parties.inp, stop.inp and xdmfWriter.inp


## mordant_testcase.sh
This files starts the testcase (copies the Boundary_Mordant.h-File to the coressponding location, asks for number of processors for run etc.)


## mordant_testcase.py
Carries out the final validation (is called automatically by the mordant_testcase.sh-script.


## reference_data.dat
Contains the data against which the data that are calculated during testrun are compared.

Important: The original testcase according to Mordant et al. requires a relativbely large domain, leading to a long computation time. After is twas verified that the Parties-Code can reproduce the experimental data from Mordant, the testcase was implified in order carry out the simulation much faster. The data contained in "reference_data.dat" are the ones from the simplified setup.


## parties_orig.inp, p_mobile_orig.inp
Represent the two "original" input files for the testcase, meaning according to the mordant testcase (not simplified).


## reference_data_orig.inp
Contains the original simulation data that was obtained by the original Mordant-Testcase (not simplified).


## plot_reference_vs_simulation
Contains an svg-graphic that displays the original experimental data (Mordant), original simulation data (Biegert) as well as the "new" Parties simulation data.
