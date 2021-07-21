# Testcases guideline
Testcases are required to keep the consistency of the code and avoid undesireable changes in the base parts of the code 
which are responsible for output results (such as numerical solvers, discretization schemes and others).

The testcases are used in case of a pull request which includes changed or added numerical functionality.
Before merging your branch with a master branch you have to perform all the testcases simulations to check 
if your additions to the code affected the results of the testcases simulations.

It is also possible that the results could be different than before. In this case you have to clearly state 
the reason and provide a strong motivation for your changes.

In case if your updated code provides an additional functionality, it is highly encouraged to create a new test case,
which would represent a simulation where your new functionality is tested. Below you can find a section which describes the exact 
procedure how a new test case should be added (section __Add a testcase__).


## Testcases structure
The directory `<$PARTIES_home>/testcases` is assigned for the testcases related files and folders.

The following folders and files are located in `<$PARTIES_home>/testcases`:
 - The folders which correspond to a certain testcase (e.g. Poiseuille_Testcase(link)).
 - The main script file which runs and checks all testcases runtestcases.sh(link).
 - The folder which contains validated results of testcases Testcase_Results(link).

The folder with a certain test case contains the following files and folders:
 - The setup files which are used to perform a certain testcase (e.g. parties.inp(link), Boundary.h(link), p_fixed.inp(link)).
 - The individual script file which is executed to run and check a corresponding testcase (e.g. Poiseuille.sh(link)).
 - Files which contain the data for initial validation of a testcase. These are usually the files related to the output data from anlitical solution or experimental results.


## Run testcases
### All testcases
To run all test cases, go to the folder testcase(add link), open terminal and run the script which includes all test cases:
```
$ cd <$PARTIES_home>/testcases/
$ ./runtestcases.sh > logfile.log
```
The script `runtestcases.sh` includes the list of all test cases and executes them in a sequential manner.
The status of each test case is written in the logfile `logfile.log`.
Therefore if the status of every test case is "passed" the code could be submitted for a pull request.
In case if some of test cases did not passed you might think about the reasons of it.
If you have a certain explanation for it, for example you implemented a new method for solving PDE, try to compare your results with analitical or experimental results.
In case if you dont have an explanation why results are different - try to double check your code and find out why the results are different.

__Note :__ the intermediate results for each test case are deleted after the passing criteria is checked. 

### Certain testcase 
If you want to analize the simulation data for a certain testcase in more details, run a script for this testcase (e.g. Poiseuille.sh(link)).
Each individual script uses the simulation setup files from corresponding folder and compiles the PARTIES code based on them.
Simulations starts after the compilation the code. Finally, the ouput data is compared with the results file from the same simulation.
The simulation results are stored in the work folder which is automatically created in `<$PARTIES_home>/testcases/ (check it)`.


## Add a testcase
A new testcase could be added using the following procedure:

1. Create a folder with your testcase into `<$PARTIES_home>/testcases/`.
2. Add all necessary files which would define the initial and boundary condition of your testcase. Use Poiseuille_Testcase and Couette_Testcase as a reference.
3. Into your testcase folder create a local script which would compile the code, run simulation and check (and print) the passing criteria (again, use Poiseuille_Testcase and Couette_Testcase as a reference)
4. Modify the `runtestcases.sh` in the following way:
- Create a variable pointing to your local script file <br> `NAME_OF_YOUR_TESTCASE="$test_home_path/NAME_OF_YOUR_TESTCASE_Testcase/NAME_OF_YOUR_TESTCASE_testcase.sh"`
- Add comands to print and execute your local script file <br>
```
printf '\n\n'
cd NAME_OF_YOUR_TESTCASE_Testcase/
. "$NAME_OF_YOUR_TESTCASE"
```
- Add comands to remove the working directory for your testcase <br> `rm -rf NAME_OF_YOUR_TESTCASE_run`

After modifications the `runtestcases.sh` file should look like this:
```
#!/bin/bash
# Shell script to run selected test flows with the PARTIES fluid solver.
# v1
# Ralph George: r.george@tu-braunschweig.de

# Current directory will be <PARTIES_home>/testcases/
echo '----------**Starting PARTIES Testcases**----------'
test_home_path=`pwd`

# Path to executables
POISEUILLE_FLOW="$test_home_path/Poiseuille_Testcase/poiseuille_flow_testcase.sh"
COUETTE_FLOW="$test_home_path/Couette_Testcase/couette_flow_testcase.sh"
NAME_OF_YOUR_TESTCASE="$test_home_path/NAME_OF_YOUR_TESTCASE_Testcase/NAME_OF_YOUR_TESTCASE_testcase.sh"
###########################################################################
#                          Testcase Execution                             #
###########################################################################

# Comment out unrequired flows
printf '\n\n'
cd Poiseuille_Testcase/
. "$POISEUILLE_FLOW"

printf '\n\n'
cd Couette_Testcase/
. "$COUETTE_FLOW"

printf '\n\n'
cd NAME_OF_YOUR_TESTCASE_Testcase/
. "$NAME_OF_YOUR_TESTCASE"

###########################################################################
#                                Clean up                                 #
###########################################################################
cd $test_home_path
rm -rf Poiseuille_run
rm -rf Couette_run
rm -rf NAME_OF_YOUR_TESTCASE_run
```
