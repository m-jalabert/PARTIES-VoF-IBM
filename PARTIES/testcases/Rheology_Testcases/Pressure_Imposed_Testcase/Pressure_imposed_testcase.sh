#!/bin/bash
# Shell script to run a test for Pressure Imposed Testcase with the PARTIES fluid solver.
# v1

# Current directory will be /testcases/Rheology_Testcase/Pressure_Imposed_Testcase
echo '----------**Start Pressure Imposed Testcase**----------'
home_path=`pwd`

###########################################################################
#                         Set no of processors                            #
###########################################################################

if [ "$1" != "" ]; then
    # no of processors is passed
    echo "Executing on $1 processors"
else
    # no of processors is not passed
    read -p 'Enter number of processors: ' nproc
    set -- "$nproc"
    echo "Executing on $1 processors"
fi

###########################################################################
#              Placing the required files in working dir                  #
###########################################################################
cd ../../
mkdir Pressure_Imposed_run
chmod 777 Pressure_Imposed_run
cp $home_path/*.* Pressure_Imposed_run/
cd Pressure_Imposed_run
work_path=`pwd`

###########################################################################
#                        Setting the Boundary.h                           #
###########################################################################
cp Boundary_PI.h ../../src/Include
cd ../../src/Include
boundary_path=`pwd`
mv Boundary.h Boundary_ORIG.h
mv Boundary_PI.h Boundary.h

###########################################################################
#                          Running parties.sh                             #
###########################################################################
cd ../../
echo 'START: Compiling executable'
printf 'make clean\nmake\n'
make clean > make_Pressure_Imposed.log
make >> make_Pressure_Imposed.log 2>&1
echo 'END: Compiling executable'

cp parties $work_path
mv make_Pressure_Imposed.log $work_path
cd $boundary_path
rm -rf Boundary.h
mv Boundary_ORIG.h Boundary.h

cd $work_path
echo 'START: Flow simulation'
mpirun -np $1 parties > output.log 2>&1
echo 'END: Flow simulation'

###########################################################################
#                        Validation 1 paticle                             #
###########################################################################
#  Run the python-script
python3 Pressure_Imposed.py
STATUS=$?
if [ $STATUS -ne 0 ]; then
        (exit $STATUS)
fi



###########################################################################
#                        Validation 2 fluid                           #
###########################################################################

###########################################################################
#                           Reading Data_*.h5                             #
###########################################################################
# Obtains the velocities at the selected plane of the grid
h5dump -d "/u" -s "45,0,45" -c "46,91,1" -w 0 Data_1.h5 > velo.dat

# File manipulation: Velocities placed in a file !!Omitted the last entry!!
header_line_no=`grep -n "\<DATA\>" velo.dat | gawk '{print $1}' FS=":"`
echo '      U(y)' > numerical_velo.dat
awk "NR==$((header_line_no+1)), NR==$((header_line_no+90))" velo.dat | gawk '{print $2}' FS=": " | gawk '{print $1}' FS="," >> numerical_velo.dat

###########################################################################
#                                L2 Norm                                  #
###########################################################################
mv numerical_velo_verified_PI.dat numerical_velo_verified.dat
g++ l2norm.c -o l2norm.sh
./l2norm.sh

###########################################################################
#                                Post proc                                #
###########################################################################
# Replace file in Results folder if passed
if grep -Fxq PASS absolute_error.dat
then
    # "Test Passed"
    # delete the old verified numerical soln
    # rename the new numerical soln as the verified soln
    rm -rf numerical_velo_verified.dat
    mv numerical_velo.dat numerical_velo_verified_PF.dat
    mv absolute_error.dat absolute_error_PF.dat
    # Clean up
    rm -rf velo.dat
    cd $home_path/../../
    (exit 0)
    # cp numerical_velo_verified_PF.dat l2norm_PF.dat ../Testcase_Results/
else
    # "Test Failed"
    # old verified numerical soln, remains as is
    # new numerical soln remains as is
    mv numerical_velo_verified.dat numerical_velo_verified_PF.dat
    mv numerical_velo.dat numerical_velo_failed_PF.dat
    mv absolute_error.dat absolute_error_failed_PF.dat
    # Clean up
    rm -rf velo.dat
    cd $home_path/../../
    (exit 1)
fi
