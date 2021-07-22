# xdmf_Writer

The xdmf_Writer creates Reader_*.xmf files to visualize the simulation data, e.g. in ParaView or Visit. The [xdmf_Writer directory](https://github.com/metialex/PARTIES/blob/kleischmann_xdmfWriter/PARTIES/matlab/xdmf_Writer) contains the main script <em>xdmf_Writer.m</em> as well as the functional scripts <em>write_Reader_c_xmf.m</em>, <em>write_Reader_Particle_xmf.m</em>, <em>write_Reader_Scalar_Velocity_xmf.m</em>, and <em>write_Reader_Vector_Velocity_xmf.m</em>. <br>

The following specification deals with the main script <em>xdmf_Writer.m</em>.

<b> Prerequisite: </b> <em>Data_*.h5</em> files; the simulation must have already run. <br> 
<b> Input: </b> Enter the path to the simulation file in which the <em>Data_*.h5</em> files are located. <br>
<pre>
%% Input
% enter path to simulation folder
% if the whole file "xdmf_Writer" has been copied to the simulation folder,
% set: path = '..';
path = '/home/Desktop/PARTIES'; 
</pre>

Moreover, decide if you want to consider all existing h5-files, means from the first until the last time step, or only a specific range. If a specific is wanted, provide the start- and end-h5 file.
<pre>
% Do you want to consider all h5-files (timesteps)?
% 1 = yes, 0 = no
all_files = 1; 

% If no, please provide start- and end-file:
start_file = 0;
end_file = 5;
</pre>

<b> Process: </b> The Reader_*.xmf are created by applying the rspective functional scripts.

<b> Output: </b> Reader_c.xmf, Reader_u.xmf, Reader_v.xmf, and Reader_w.xmf are created in any case. Reader_p_mobile.xmf and/or Reader_p_fixed.xmf are created in case that a respective particle is specified. Reader_vector.xmf is created in case the Vector_*.h5 files exist. <br>
