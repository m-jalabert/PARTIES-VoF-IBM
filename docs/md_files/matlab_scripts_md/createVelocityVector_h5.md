# createVelocityVector_h5.m

The present matlab script creates cell centered velcoity vectors in each domain cell based on the scalar velocities u, v, and w. <br>

<b> Prerequisite: </b> <em>Data_*.h5</em> files; the simulation must have already run. <br> 
<b> Input: </b> Enter the path to the simulation file in which the <em>Data_*.h5</em> files are located. <br>
<pre>
%% Input
% enter path to simulation folder
% if matlab script has been copied to simulation folder, set: path = '';
path = '/home/Desktop/PARTIES';
</pre>

<b> Process: </b> The face velocities of each cell are averaged in each coordinate direction to calculate the cell centered velocities. <br>
<pre>
for i = 1 : Nx_i
  for j = 1 : Ny_i
    for k = 1 : Nz_i
      % calculate the centered velocities for u, v, and w
      velocity_matrix_4D(1, i, j, k) = 0.5 * (double(u(i, j, k)) + double(u(i+1, j, k)));
      velocity_matrix_4D(2, i, j, k) = 0.5 * (double(v(i, j, k)) + double(v(i, j+1, k)));
      velocity_matrix_4D(3, i, j, k) = 0.5 * (double(w(i, j, k)) + double(w(i, j, k+1)));
    end %k
  end %j
end %i
</pre>

<b> Output: </b> The velocity vectors are saved as <em>Vector_*.h5</em> files. Each <em>Vector_*.h5</em> file contains the timestep specific velocity vectors and the group "grid", which includes the total number of grid points in each direction (NX, NY, NZ) as well as the cell centered coordinates (xc, yc, zc). <br>





