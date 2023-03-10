import math
import os
import sys
import csv

### Load reference data ###
f           = open(os.path.join(sys.path[0], "reference_data.dat"), 'r')
csv_f       = csv.reader(f, delimiter=',')
reference   = []
for row in csv_f:
     reference.append(row[1:2])                                 # Get the third column (x-position)
f.close()
flatList    = [item for elem in reference for item in elem]     # Flatten the nested list obtained by for-loop
reference   = [float(i) for i in flatList]                      # Convert each list element into a float

### Load simulation data ###
f           = open(os.path.join(sys.path[0], "mobile.dat"), 'r')
csv_f       = csv.reader(f, delimiter=',')
simulation   = []
for row in csv_f:
     simulation.append(row[1:2])                                # Get the third column (x-position)
f.close()
flatList    = [item for elem in simulation for item in elem]    # Flatten the nested list obtained by for-loop
simulation  = [float(i) for i in flatList]                      # Convert each list element into a float


### Create delta list of reference and simulation data ###
delta = []
for i in range(0, len(reference)):                              # Create list of delta between reference and simulation values
    if reference[i] != 0.0:
        delta.append(abs(reference[i]-simulation[i])/abs(reference[i]))

### Check for test condition ###
result = all(ele < 1e-12 for ele in delta)                        # Check whether each element of delta list is smaller than 1e-12
if result == True:
    print('\nTest passed')
    print('----------**End Volume Imposed Testcase part I**----------\n')
    sys.exit(0)
else:
    print('\nTest failed.')
    print('\n----------**End Volume Imposed Testcase**----------')
    sys.exit(1)
