#!/bin/bash
# Shell script to run a test for Poiseuilli Flow with the PARTIES fluid solver.
# v1
# Ralph George: r.george@tu-braunschweig.de

# Current directory will be testcases/<Flow Type>
echo 'Starting PARTIES test run: Poiseuilli flow'
echo 'Version 1.0'

home_path=`pwd`
echo 'Home: '$home_path

count=0

###########################################################################
#             Selecting the correct boundary condition                    #
###########################################################################
# cd ../Compilation/src/Include/
cd ../../src/Include/
cp Boundary.h{,_BKP}

# Change all the current BCs to undef state
sed -i '4,$ s/\<define\>/undef/' Boundary.h

# Implement the correct BC for Poiseuilli Flow
for i in `grep -n PERIODIC_NOSLIP_BOX Boundary.h | gawk '{print $1}' FS=":"`
do
    let "count+=1"
    if [ $count -eq 1 ]
    then
        sed -i "$i s/undef/define/" Boundary.h
    elif [ $count -eq 2 ]
    then
        sed -i "$((i+1)) s/undef/define/" Boundary.h
        sed -i "$((i+2)) s/undef/define/" Boundary.h
        sed -i "$((i+3)) s/undef/define/" Boundary.h
        sed -i "$((i+4)) s/undef/define/" Boundary.h
    fi
done

###########################################################################
#               Enabling CONSTANT_MASSFLUX and CG_SOLVE                   # 
###########################################################################
cm_line_no=`grep -n CONSTANT_MASSFLUX Boundary.h | gawk '{print $1}' FS=":"`
sed -i "$cm_line_no s/undef/define/" Boundary.h

cg_line_no=`grep -nw CG_SOLVE Boundary.h | gawk '{print $1}' FS=":" | head -1`
sed -i "$cg_line_no s/undef/define/" Boundary.h

###########################################################################
#                          Running parties.sh                             #
###########################################################################
# All input files(.inp) with the required values are present
cd ../../
make clean
make

cp parties $home_path
mpirun -np 2 parties > output.log

###########################################################################
#                          Analytical solution                            #
###########################################################################
# Read dat file from analytical solution
# analytical soluntion already present for given flow

###########################################################################
#                           Reading Data_5.h5                             #
###########################################################################
# Obtains the velocities at the selected plane of the grid
h5dump -d "/u" -s "0,0,150" -c "1,101,1" -w 0 Data_10.h5 > velo.dat

# File manipulation: Velocities placed in a file !!Omitted the last entry!!
header_line_no=`grep -n "\<DATA\>" velo.dat | gawk '{print $1}' FS=":"`
echo '      U(y)' > numerical_velo_PF.dat
awk "NR==$((header_line_no+1)), NR==$((header_line_no+100))" velo.dat | gawk '{print $2}' FS=": " | gawk '{print $1}' FS="," >> numerical_velo_PF.dat

###########################################################################
#                              Comparision                                #
###########################################################################
# Gives numerical minus analytical: Absolute value
awk 'FNR==NR { _a[FNR]=$1;} NR!=FNR { $1 -= _a[FNR]; {sub("^-", "", $1); if ($1 > 0.02) print $1,"FNR="FNR};  }'  analytical_velo_PF.dat numerical_velo_PF.dat


###########################################################################
#                                Clean up                                 #
###########################################################################
rm -rf velo.dat