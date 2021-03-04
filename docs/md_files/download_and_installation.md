# Downloading and Unpacking the Code

[PARTIES_Install](https://github.com/metialex/PARTIES/tree/master/PARTIES_Install) is included in the repository and contains the required libraries (<em>OpenMPI</em>, <em>HDF5</em>, and <em>FFTW</em>) for installing or compiling, respectively, the PARTIES code. Therefore, you should firstly clone the repository to your local machine.

<pre>
git clone <em>URL</em>
</pre>

The libraries are in <em>.tar</em> and <em>.tar.gz</em> format and need to be unpacked. Unpack the files to a local directory, which we will call <em>your_directory</em>, and create an additional file, e.g. <em>PARTIES_Libs</em>, in the same directory. Please be aware that you have to insert the path of your directory for <em>your_directory</em> in all following commands.

<pre>
tar -xf <em>*.tar*</em> -C <em>your_directory</em>
cd <em>your_directory</em>
mkdir <em>PARTIES_Libs</em>
</pre>

# Installation

The three libraries <em>OpenMPI</em>, <em>HDF5</em>, and <em>FFTW</em> have to be installed successively.

## OpenMPI

Go into the unpacked <em>OpenMPI</em>-file within <em>your_directory</em> and run the following commands one after the other. 

<pre>
./configure --prefix=$HOME/<em>your_directory</em>/PARTIES_Libs/
make all install
</pre>

After the installation, you have to add the BUILTIN command "export..." to your bash script, which will refer to the previously installed /bin-file of <em>OpenMPI</em>. Therefore, open the bash script (e.g. with the text editor <em>nano</em>), add the cited line, save and close the document, and source the bash script as shown in the following commands:

<pre>
nano ~/.bashrc
<em>add</em> the line:   export PATH=$HOME/<em>your_directory</em>/PARTIES_Libs/bin:$PATH
<em>save</em> and close
source ~/.bashrc
</pre>

## HDF5

Change into the unpacked <em>HDF5</em>-file and run the commands below. Please double check that you have replaced all <em>your_directory</em> references with  your path.


