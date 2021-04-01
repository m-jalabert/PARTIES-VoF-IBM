#!/bin/bash
# Shell script to run a test for Poiseuille Flow with the PARTIES fluid solver.
# v1
# Ralph George: r.george@tu-braunschweig.de

# Current directory will be <PARTIES_home>/testcases/Poiseuilli_Testcase/
echo '----------**Start Mordant  Testcase**----------'
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
cd ../
mkdir Mordant_run
chmod 777 Mordant_run
cp $home_path/*.* Mordant_run/
cd Mordant_run
work_path=`pwd`

###########################################################################
#                        Setting the Boundary.h                           #
###########################################################################
cp Boundary_Mordant.h ../../src/Include
cd ../../src/Include
boundary_path=`pwd`
mv Boundary.h Boundary_ORIG.h
mv Boundary_Mordant.h Boundary.h

###########################################################################
#                          Running parties.sh                             #
###########################################################################
cd ../../
echo 'START: Compiling executable'
printf 'make clean\nmake\n'
make clean > make_mordant.log
make >> make_mordant.log
echo 'END: Compiling executable'

cp parties $work_path
mv make_mordant.log $work_path
cd $boundary_path
rm -rf Boundary.h
mv Boundary_ORIG.h Boundary.h

cd $work_path
echo 'START: Flow simulation'
mpirun -np $1 parties > output.log
echo 'END: Flow simulation'

###########################################################################
#                           Validation                             #
###########################################################################
#  Run the python-script
python mordant_testcase.py

