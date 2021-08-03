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

###########################################################################
#                          Checking mobile.dat                            #
###########################################################################
echo 'START: Data comparison'
python3 gondret_st_152_compare.py
echo 'END: Data comparison'
