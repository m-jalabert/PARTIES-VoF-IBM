#!/bin/bash
#####################################################################
#
# Project 	:	PARTIES flow solver
# Program name	:	install_parties.sh
# Author	:	Thomas Köllner
# Date created	:	January 30th 2019
# Purpose	: 	Script to install PARTIES. No special rights of the user account are 				necessary. Onyl a C-compiler has to be present; I assume the GNU compiler 				collection can be found in the $PATH and that the sources zip §PARTIES_ALL 				can be found. 
# Revision	:	
# Date		Author	Ref	Revision (Date in DDMMYYYY)
# 01022019	TK	1	Initial Creation; dowloading the files from our Gdrive might be unstable and prone to errors due to changes in googles website. This step might be done manually or the files shoul be hostes on our own server. 
###################################################################### 

BASE_DIR=${HOME}/Software/
INSTALL_DIR=${BASE_DIR}PARTIES_Libs/  ## this can be adapted to your environment

FFTW_SOURCE=${BASE_DIR}fftw-3.3.8
HDF5_SOURCE=${BASE_DIR}hdf5-1.8.20
MPI_SOURCE=${BASE_DIR}openmpi-3.1.3

MPI_BUILD=${INSTALL_DIR}openmpi-3.1.3

SWITCH_BUILD_HDF5=1
SWITCH_BUILD_MPI=1
SWITCH_BUILD_FFTW=1

#################

echo -n "This scripts builds the libraries for PARTIES. All files in ${INSTALL_DIR} will be deleted or if the directory is not present, it will be created. Do you wish to continue (y/n)? "
read answer


if [ "$answer" != "${answer#[Yy]}" ]; then
	echo "YES \n"

	if [ ! -d "${BASE_DIR}" ]; then
		mkdir ${BASE_DIR}
	fi

	if [ ! -d "${INSTALL_DIR}" ]; then
	  	mkdir $INSTALL_DIR
	else
		mv   "$INSTALL_DIR"  "${INSTALL_DIR%/}_old" 
		echo "The old installation was moved to ${INSTALL_DIR%/}_old \n" 

		mkdir $INSTALL_DIR
	fi


else
    	echo No
    	exit 0
fi

####### Now we download the source files and unzip them

cd ${BASE_DIR}

if [ -f ${HDF5_SOURCE}.tar  ]
then
	tar -xf ${HDF5_SOURCE}.tar
else
echo "\n Download the HDF5 library from our Google Drive \n" ##  if this does not work automatically the files have to be downloaded manually
	curl -L -o ${HDF5_SOURCE}.tar "https://drive.google.com/uc?export=download&id=1xrz1N96A3Kx7xymepj6B8-SpGGUb7Le2" # download 
	tar -xf ${HDF5_SOURCE}.tar  
fi



if [ -f ${FFTW_SOURCE}.tar.gz  ]
then
	tar -xf ${FFTW_SOURCE}.tar.gz
else
echo "\nDownload the FFTW library from our Google Drive \n" ##  if this does not work automatically the files have to be downloaded manually
	curl -L -o ${FFTW_SOURCE}.tar.gz "https://drive.google.com/uc?export=download&id=117zEaDE6c7djJv5Vw6OJ7hDerqSCmdhw" # download 
	tar -xf ${FFTW_SOURCE}.tar.gz
fi



if [ -f ${MPI_SOURCE}.tar.gz  ]
then
	tar -xf ${MPI_SOURCE}.tar.gz
else
echo "\n Download the MPI library from our Google Drive \n" ##  if this does not work automatically the files have to be downloaded manually

	filename="${MPI_SOURCE}.tar.gz"
	fileid="1Mv3DcDLh79TFosIYwAg0_3FlY67IOPWf"
	
curl -c ./cookie -s -L "https://drive.google.com/uc?export=download&id=${fileid}" > /dev/null
curl -Lb ./cookie "https://drive.google.com/uc?export=download&confirm=`awk '/download/ {print $NF}' ./cookie`&id=${fileid}" -o ${filename} 



	#wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=$FILEID' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=$FILEID" -O $FILENAME #&& rm -rf /tmp/cookies.txtman curl
## this is the save way to download things see https://gist.github.com/tanaikech
#	query=′curl -c ./cookie.txt -s -L "https://drive.google.com/uc?export=download&id=${fileid}" | pup 'a#uc-download-link attr{href}' | sed -e 's/amp;//g'′
#	curl -b ./cookie.txt -L -o ${filename} "https://drive.google.com${query}"
#	curl -L -o ${MPI_SOURCE}.tar.gz "https://drive.google.com/uc?export=download&id=${fileid}" # download 

	tar -xf ${MPI_SOURCE}.tar.gz
fi

##### end downloading and unzipping


#read -p "\n Please, press [enter] to build the OpenMPI library \n"
#### Building OpenMPI
if [ ${SWITCH_BUILD_MPI} -eq "1" ] 
then
echo "\n............................ Okay building the MPI library.....................................\n"
cd ${MPI_SOURCE}
./configure --prefix=${INSTALL_DIR}
make all install
fi
## for Mac it must be rpath,$  so change equal sign to comma
#### Building HDF5
if [ ${SWITCH_BUILD_HDF5} -eq "1" ] 
then
echo "\n............................ Okay building the HDF5 library.....................................\n"
cd ${HDF5_SOURCE}
./configure --prefix=${INSTALL_DIR} --enable-parallel CC="${INSTALL_DIR}bin/mpicc" CXX="${INSTALL_DIR}bin/mpicxx" LDFLAGS="-L${INSTALL_DIR}lib -Wl,-rpath=${INSTALL_DIR}lib -fPIC" CPPFLAGS="-I${INSTALL_DIR}include" CFLAGS="-I${INSTALL_DIR}include"
make 
make check
make install
make check-install
fi

#### Building FFTW
if [ ${SWITCH_BUILD_FFTW} -eq "1" ] 
then
echo "\n............................ Okay building the FFTW library.....................................\n"
cd ${FFTW_SOURCE}

./configure --prefix=${INSTALL_DIR} --enable-mpi CFLAGS="-O3" MPICC="${INSTALL_DIR}bin/mpicc" LDFLAGS="-L${INSTALL_DIR}lib/openmpi -Wl,-rpath=${INSTALL_DIR}lib/openmpi -fPIC" CPPFLAGS="-I ${INSTALL_DIR}include"
make
make install
make installcheck
fi


