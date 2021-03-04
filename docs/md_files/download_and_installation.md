# Downloading and Unpacking the Code

[PARTIES_Install](https://github.com/metialex/PARTIES/tree/master/PARTIES_Install) is included in the repository and contains the required libraries (<em>OpenMPI</em>, <em>HDF5</em>, and <em>FFTW</em>) for installing or compiling, respectively, the PARTIES code. Therefore, you should firstly clone the repository to your local machine.

<pre>
git clone <em>URL</em>
</pre>

The libraries are in <em>.tar</em> and <em>.tar.gz</em> format and need to be unpacked. Unpack the files to a local directory, which we will call <b>your_directory</b>, and create an additional file, e.g. <em>PARTIES_Libs</em>, in the same directory. Please note that you have to specify the path of your directory in all following commands for <b>your_directory</b>.

<pre>
tar -xf <em>*.tar*</em> -C <b>your_directory</b>
cd <b>your_directory</b>
mkdir <em>PARTIES_Libs</em>
</pre>

# Installation

The three libraries <em>OpenMPI</em>, <em>HDF5</em>, and <em>FFTW</em> have to be installed successively.

## OpenMPI

Go into the unpacked <em>OpenMPI</em> file within <b>your_directory</b> and run the following commands one after the other. 

<pre>
./configure --prefix=$HOME/<b>your_directory</b>/PARTIES_Libs/
make all install
</pre>

After the installation, you have to add the BUILTIN command "export..." to your bash script, which will refer to the previously installed /bin-file of <em>OpenMPI</em>. Therefore, open the bash script (e.g. with the text editor <em>nano</em>), add the cited line, save and close the document, and source the bash script as shown in the following commands:

<pre>
nano ~/.bashrc
<em>add the line:</em>   export PATH=$HOME/<b>your_directory</b>/PARTIES_Libs/bin:$PATH
<em>save and close</em>
source ~/.bashrc
</pre>

## HDF5

Change into the unpacked <em>HDF5</em> file and run the commands below. Please double check that you have specified all <b>your_directory</b> references with your path.

<pre>
./configure --prefix=$HOME/<b>your_directory</b>/PARTIES_Libs/ --enable-parallel CC="$HOME/<b>your_directory</b>/PARTIES_Libs/bin/mpicc" CXX="$HOME/<b>your_directory</b>/PARTIES_Libs/bin/mpicxx" LDFLAGS="-L$HOME/<b>your_directory</b>/PARTIES_Libs/lib -fPIC" CPPFLAGS="-I$HOME/<b>your_directory</b>/PARTIES_Libs/include"
make
make check
make install
make check-install
</pre>

## FFTW

Change into the unpacked <em>FFTW</em> file and run the commands below. Similar to the installation of <em>HDF5</em>, please double check that you have specified all <b>your_directory</b> references with your path.

<pre>
./configure --prefix=$HOME/<b>your_directory</b>/PARTIES_Libs/ --enable-mpi CFLAGS="-O3" MPICC="$HOME/<b>your_directory</b>/PARTIES_Libs/bin/mpicc" LDFLAGS="-L$HOME/<b>your_directory</b>/PARTIES_Libs/lib/openmpi -Wl,-rpath=$HOME/<b>your_directory</b>/PARTIES_Libs/lib/openmpi -fPIC" CPPFLAGS="-I $HOME/<b>your_directory</b>/PARTIES_Libs/include"
make
make install
make installcheck
</pre>


# Compilation

The installation is now completed, so you can compile the PARTIES code as an example to check if it works. To do this, please go to the <em>PARTIES</em> file in your local Github repository and execute the following:

<pre>
make clean
make
</pre>

In case of a successful compilation, you are now able to run a simulation as shown in [Simulation setup](https://github.com/metialex/PARTIES/blob/master/docs/md_files/simulation_setup.md).
