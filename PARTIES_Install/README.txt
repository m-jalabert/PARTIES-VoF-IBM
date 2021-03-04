1. 	make the shell script executeable,if necessary				
	=>   'chmod +x install_script.sh'

2.	run the installation script
	=>   './install_script.sh'

3.	add the following line into ~/.bashrc
	=>   export PATH=$HOME/Software/PARTIES_Libs/bin:$PATH

4.	source bashrc
	=>   'source ~/.bashrc'

5.	run the following commands:
	=>   cd ~/Software/hdf5-1.8.20/
	=>   ./configure --prefix=$HOME/Software/PARTIES_Libs/ --enable-parallel CC="$HOME/Software/PARTIES_Libs/bin/mpicc" CXX="$HOME/Software/PARTIES_Libs/bin/mpicxx" LDFLAGS="-L$HOME/Software/PARTIES_Libs/lib -fPIC" CPPFLAGS="-I$HOME/Software/PARTIES_Libs/include"
	=>   make
	=>   make check
	=>   make install
	=>   make check-install

	=>   cd ~/fftw-3.3.8/
	=>   ./configure --prefix=$HOME/Software/PARTIES_Libs/ --enable-mpi CFLAGS="-O3" MPICC="$HOME/Software/PARTIES_Libs/bin/mpicc" LDFLAGS="-L$HOME/Software/PARTIES_Libs/lib/openmpi -Wl,-rpath=$HOME/Software/PARTIES_Libs/lib/openmpi -fPIC" CPPFLAGS="-I $HOME/Software/PARTIES_Libs/include"
	=>   make
	=>   make install
	=>   make installcheck

6.	go into $HOME/Software/PARTIES_new/ and run the following commands:
	=>   'make clean'
	=>   'make'
