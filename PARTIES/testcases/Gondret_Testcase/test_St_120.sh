#!/bin/bash
# Shell script to run a test for particle contact (Gondret test case) with the PARTIES fluid solver.
# v1
# Alexander Metelkin: a.metelkin@tu-braunschweig.de

# Current directory will be <PARTIES_home>/testcases/Gondret_Testcase/
echo '----------**Collision test case St = 120**----------'
cd test_St_120
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
mkdir test_St_120_run
chmod 777 test_St_120_run
cp $home_path/*.* test_St_120_run/
cd test_St_120_run
work_path=`pwd`

###########################################################################
#                        Setting the Boundary.h                           #
###########################################################################
cp Boundary_St_120.h ../../src/Include
cd ../../src/Include
boundary_path=`pwd`
mv Boundary.h Boundary_ORIG.h
mv Boundary_St_120.h Boundary.h


###########################################################################
#                          Running parties.sh                             #
###########################################################################
cd ../../
echo 'START: Compiling executable'
printf 'make clean\nmake\n'
make clean > make_St120.log
make >> make_St120.log
echo 'END: Compiling executable'

cp parties $work_path
mv make_St120.log $work_path
cd $boundary_path
rm -rf Boundary.h
mv Boundary_ORIG.h Boundary.h

cd $work_path
echo 'START: Flow simulation'
mpirun -np $1 parties > output.log
echo 'END: Flow simulation'

###########################################################################
#                          Checking mobile.dat                            #
###########################################################################
echo 'START: Data comparison'
python3 test_st_120_compare.py
STATUS=$?
cd $home_path/../..
echo 'END: Data comparison'
(exit $STATUS)