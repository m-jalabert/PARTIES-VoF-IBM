import math
import os
import sys
import csv
from termcolor import colored

### Load reference data ###
f           = open(os.path.join(sys.path[0], "../Gondret_Testcase/Reference_data/mobile_st_27.dat"), 'r')
csv_f       = csv.reader(f, delimiter=',')
reference   = []
for row in csv_f:
     reference.append(row[3:4])                                 # Get the fourth column
f.close()
flatList    = [item for elem in reference for item in elem]     # Flatten the nested list obtained by for-loop
reference   = [float(i) for i in flatList]                      # Convert each list element into a float

### Load simulation data ###
f           = open(os.path.join(sys.path[0], "mobile.dat"), 'r')
csv_f       = csv.reader(f, delimiter=',')
simulation   = []
for row in csv_f:
     simulation.append(row[3:4])                                # Get the fourth column
f.close()
flatList    = [item for elem in simulation for item in elem]    # Flatten the nested list obtained by for-loop
simulation  = [float(i) for i in flatList]                      # Convert each list element into a float


### Create delta list of reference and simulation data ###
delta = []
for i in range(0, len(reference)):                              # Create list of delta between reference and simulation values
    delta.append(abs(reference[i]-simulation[i]))
    if (delta[i] > 1e-12):
        print(colored(('\nTest failed. Difference {} in line {}'.format(delta[i],i+1)), 'red'))
        exit()

print(colored('\nTest passed.','green'))
