# unpack all files
cd ~/Software
tar -xf ~/InstallationFiles/PARTIES_Install/openmpi-3.1.3.tar.gz
tar -xf ~/InstallationFiles/PARTIES_Install/hdf5-1.8.20.tar
tar -xf ~/InstallationFiles/PARTIES_Install/fftw-3.3.8.tar.gz
tar -xf ~/InstallationFiles/PARTIES_Install/PARTIES_new.tar.gz


mkdir PARTIES_Libs


# OpenMPI
cd openmpi-3.1.3/
./configure --prefix=$HOME/Software/PARTIES_Libs/ 
make all install 



