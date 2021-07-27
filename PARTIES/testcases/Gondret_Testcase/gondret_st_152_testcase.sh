#!/bin/bash
# Shell script to run a test for particle contact (Gondret test case) with the PARTIES fluid solver.
# v1
# Alexander Metelkin: a.metelkin@tu-braunschweig.de

# Current directory will be <PARTIES_home>/testcases/Couette_Testcase/
echo '----------**Gondret test case St = 152**----------'
cd St_152
home_path=`pwd`
cd ../
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
#             Placing the required files in working dir                   #
###########################################################################
cd ../
mkdir St_152_run
chmod 777 St_152_run
cp $home_path/*.* St_152_run/
cd St_152_run
work_path=`pwd`

###########################################################################
#                        Setting the Boundary.h                           #
###########################################################################
cp Boundary_G152.h ../../src/Include
cd ../../src/Include
boundary_path=`pwd`
mv Boundary.h Boundary_ORIG.h
mv Boundary_G152.h Boundary.h


###########################################################################
#                          Running parties.sh                             #
###########################################################################
cd ../../
echo 'START: Compiling executable'
printf 'make clean\nmake\n'
make clean > make_G152.log
make >> make_G152.log
echo 'END: Compiling executable'

cp parties $work_path
mv make_G152.log $work_path
cd $boundary_path
rm -rf Boundary.h
mv Boundary_ORIG.h Boundary.h

cd $work_path
echo 'START: Flow simulation'
mpirun_o -np $1 parties > output.log
echo 'END: Flow simulation'
:'
###########################################################################
#                           Reading Data_*.h5                             #
###########################################################################
# Obtains the velocities at the selected plane of the grid
h5dump -d "/u" -s "0,0,150" -c "1,101,1" -w 0 Data_1.h5 > velo.dat

# File manipulation: Velocities placed in a file !!Omitted the last entry!!
header_line_no=`grep -n "\<DATA\>" velo.dat | gawk '{print $1}' FS=":"`
echo '      U(y)' > numerical_velo.dat
awk "NR==$((header_line_no+1)), NR==$((header_line_no+101))" velo.dat | gawk '{print $2}' FS=": " | gawk '{print $1}' FS="," >> numerical_velo.dat

###########################################################################
#                                L2 Norm                                  #
###########################################################################
mv numerical_velo_verified_CF.dat numerical_velo_verified.dat
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
    mv numerical_velo.dat numerical_velo_verified_CF.dat
    mv absolute_error.dat absolute_error_CF.dat
    # cp numerical_velo_verified_CF.dat l2norm_CF.dat ../Testcase_Results/
else
    # "Test Failed"
    # old verified numerical soln, remains as is
    # new numerical soln remains as is
    mv numerical_velo_verified.dat numerical_velo_verified_CF.dat
    mv numerical_velo.dat numerical_velo_failed_CF.dat
    mv absolute_error.dat absolute_error_failed_CF.dat
fi

# Clean up
rm -rf velo.dat
cd $home_path/..
'
