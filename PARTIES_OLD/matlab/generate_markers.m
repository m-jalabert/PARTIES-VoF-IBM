function [] = generate_markers(hdf5_idx, N_pts, Nproc)
% GENERATE_MARKERS   Generates grid for visualizing sphere rotation
%    GENERATE_MARKERS(hdf5_idx, N_pts, Nproc) Creates Marker_*.h5 files
%    where * is a vector of integers specified by the input 'hdf5_idx'.
%    Within the files are contained information to plot a grid on particle
%    surfaces with a texture, as interpreted by the xdmfWriter.  'N_pts'
%    will create a sphere grid of resolution 'N_pts x N_pts'.  'Nproc' is
%    the number of processors to run on, to speed up the work.
%
%    Example:
%       generate_markers(0:50, 20, 8)
%
%    will generate files Marker_0.h5 through Marker_50.h5 with a sphere
%    mesh of 20x20 per particle, and will use 8 processors.

% Debugging to plot set of points
% [X_black, X_red, topology_black, topology_red] = generate_sphere_pts(N_pts);
% 
% figure(1), clf, hold on
% for i=1:size(topology_black,1)
% 	for j=1:4
% 		pt1 = topology_black(i,j) + 1;
% 		pt2 = topology_black(i,mod(j,4)+1) + 1;
% 		plot3([X_black(pt1,1) X_black(pt2,1)], [X_black(pt1,2) X_black(pt2,2)], ...
% 			[X_black(pt1,3) X_black(pt2,3)], 'k');
% 	end
% end
% for i=1:size(topology_red,1)
% 	for j=1:4
% 		pt1 = topology_red(i,j) + 1;
% 		pt2 = topology_red(i,mod(j,4)+1) + 1;
% 		plot3([X_red(pt1,1) X_red(pt2,1)], [X_red(pt1,2) X_red(pt2,2)], ...
% 			[X_red(pt1,3) X_red(pt2,3)], 'r--');
% 	end
% end
% 
% return

if (Nproc > 1)
	poolobj = parpool(Nproc);
	parfor count = 1:length(hdf5_idx)
		generate_markers_serial(hdf5_idx(count), N_pts)
	end
	delete(poolobj)
else
	for count = 1:length(hdf5_idx)
		generate_markers_serial(hdf5_idx(count), N_pts)
	end
end


% ******************************************************************************
% Generate markers for one HDF5 file on one processor
% ******************************************************************************
function [] = generate_markers_serial(hdf5_idx, N_pts)

% Generate set of points for a generic sphere
[X_black, X_red, topology_black, topology_red] = generate_sphere_pts(N_pts);

% Number of point and quadrilaterals
N_pts_black = size(X_black,1);
N_pts_red = size(X_red,1);
N_quad_black = size(topology_black,1);
N_quad_red = size(topology_red,1);

filename = sprintf('%s/Particle_%d.h5', pwd, hdf5_idx);
info = h5info(filename,'/');

% Find number of fixed and mobile particles
Npf = 0;
Npm = 0;
for j=1:length(info.Groups)
	group = info.Groups(j).Name;
	
	if (strcmp(group, '/mobile') || strcmp(group, '/fixed'))		
		dsinfo = h5info(filename, [group '/R']);
		if strcmp(group, '/mobile')
			Npm = dsinfo.Dataspace.Size(2);
		elseif strcmp(group, '/fixed')
			Npf = dsinfo.Dataspace.Size(2);
		end
	end
end

% Allocate space for marker points, attributes, and topologies
X_markers_b = NaN((Npf+Npm)*N_pts_black,3);
X_markers_r = NaN((Npf+Npm)*N_pts_red,3);
att_markers_b = NaN((Npf+Npm)*N_pts_black,1);
att_markers_r = NaN((Npf+Npm)*N_pts_red,1);
top_markers_b = NaN((Npf+Npm)*N_quad_black,4);
top_markers_r = NaN((Npf+Npm)*N_quad_red,4);

ID_X_b = 1:N_pts_black;
ID_X_r = 1:N_pts_red;
ID_top_b = 1:N_quad_black;
ID_top_r = 1:N_quad_red;
for j=1:length(info.Groups)

	group = info.Groups(j).Name;
	if (strcmp(group, '/mobile') || strcmp(group, '/fixed'))

		R = h5read(filename, [group '/R'])';
		X = h5read(filename, [group '/X'])';
		U = h5read(filename, [group '/U'])';
		Rotn = h5read(filename, [group '/Rotn'])';
		
		% For each particle
		for i=1:length(R)
			
			% Add marker positions
			X_markers_b(ID_X_b,:) = transform_markers(X_black, R(i), X(i,:), Rotn(i,:));
			X_markers_r(ID_X_r,:) = transform_markers(X_red, R(i), X(i,:), Rotn(i,:));
			
			% Add attributes
			att_markers_b(ID_X_b) = 0;
			att_markers_r(ID_X_r) = norm(U(i,:));
			
			% Add topology
			top_markers_b(ID_top_b,:) = topology_black + ID_X_b(1) - 1;
			top_markers_r(ID_top_r,:) = topology_red + ID_X_r(1) - 1;
			
			% Increment IDs
			ID_X_b = ID_X_b + N_pts_black;
			ID_X_r = ID_X_r + N_pts_red;
			ID_top_b = ID_top_b + N_quad_black;
			ID_top_r = ID_top_r + N_quad_red;
		end
	end
end

% Write results to HDF5 file
marker_file = sprintf('%s/Marker_%d.h5', pwd, hdf5_idx);
fprintf('Editing %s\n', marker_file)
black_group = '/marker_black';
red_group = '/marker_red';
attribute_name = '/att';
h5create(marker_file, [black_group '/X'], size(X_markers_b'));
h5write(marker_file, [black_group '/X'], X_markers_b');
h5create(marker_file, [red_group '/X'], size(X_markers_r'));
h5write(marker_file, [red_group '/X'], X_markers_r');
h5create(marker_file, [black_group attribute_name], size(att_markers_b'));
h5write(marker_file, [black_group attribute_name], att_markers_b');
h5create(marker_file, [red_group attribute_name], size(att_markers_r'));
h5write(marker_file, [red_group attribute_name], att_markers_r');
h5create(marker_file, [black_group '/topology'], size(top_markers_b'));
h5write(marker_file, [black_group '/topology'], top_markers_b');
h5create(marker_file, [red_group '/topology'], size(top_markers_r'));
h5write(marker_file, [red_group '/topology'], top_markers_r');
% h5create(marker_file, [black_group '/topology'], size(top_markers_b'), 'Datatype', 'int32');
% h5write(marker_file, [black_group '/topology'], int32(top_markers_b'));
% h5create(marker_file, [red_group '/topology'], size(top_markers_r'), 'Datatype', 'int32');
% h5write(marker_file, [red_group '/topology'], int32(top_markers_r'));



% ******************************************************************************
% Transform a set of points 'X_pts' to follow a sphere of radius 'R',
% position 'X', and rotation matrix 'Rotn'.
% ******************************************************************************
function X_transform = transform_markers(X_pts, R, X, Rotn)

Rotn = [Rotn(1:3); Rotn(4:6); Rotn(7:9)];
X_transform = Rotn * X_pts';
X_transform = bsxfun(@plus, R * X_transform', X);



% ******************************************************************************
% Generate a set of points and their topologies for a checkerboard pattern
% ******************************************************************************
function [X_black, X_red, topology_black, topology_red] = generate_sphere_pts(N_pts)

M_pts = N_pts;

[Xs,Zs,Ys] = sphere(N_pts);

% Indices to split points into patches of red and black
idx_theta = round([0, 0.25*N_pts, 0.5*N_pts, 0.75*N_pts, N_pts]) + 1;
idx_phi = round([0, 0.5*M_pts, M_pts]) + 1;

% Count number of points
N_pt_black = 0;
N_pt_red = 0;
N_quad_black = 0;
N_quad_red = 0;
for j=1:length(idx_phi)-1
	for i=1:length(idx_theta)-1
		% Number of points in phi and theta directions of subgroup
		Msub = idx_phi(j+1) - idx_phi(j) + 1;
		Nsub = idx_theta(i+1) - idx_theta(i) + 1;
		
		if mod(i+j,2) == 0
			N_pt_black = N_pt_black + Msub * Nsub;
			N_quad_black = N_quad_black + (Msub-1) * (Nsub-1);
		else
			N_pt_red = N_pt_red + Msub * Nsub;
			N_quad_red = N_quad_red + (Msub-1) * (Nsub-1);
		end
	end
end

% Allocate space
X_black = NaN(N_pt_black, 3);
X_red = NaN(N_pt_red, 3);
topology_black = NaN(N_quad_black, 4);
topology_red = NaN(N_quad_red, 4);

N_pt_black = 0;
N_pt_red = 0;
N_quad_black = 0;
N_quad_red = 0;
for j=1:length(idx_phi)-1
	for i=1:length(idx_theta)-1
		
		% Determine if we are dealing with a red or black patch
		if mod(i+j,2) == 0
			black = 1;
		else
			black = 0;
		end
		
		% Number of points in phi and theta directions of subgroup
		Msub = idx_phi(j+1) - idx_phi(j) + 1;
		Nsub = idx_theta(i+1) - idx_theta(i) + 1;
		
		if black
			ID_pt_start = N_pt_black;
			N_pt_black = N_pt_black + Msub * Nsub;
			ID_quad_start = N_quad_black;
			N_quad_black = N_quad_black + (Msub-1) * (Nsub-1);
		else
			ID_pt_start = N_pt_red;
			N_pt_red = N_pt_red + Msub * Nsub;
			ID_quad_start = N_quad_red;
			N_quad_red = N_quad_red + (Msub-1) * (Nsub-1);
		end
		
		% Coordinates of points in subgroup
		Xsub = Xs(idx_phi(j):idx_phi(j+1),idx_theta(i):idx_theta(i+1))';
		Ysub = Ys(idx_phi(j):idx_phi(j+1),idx_theta(i):idx_theta(i+1))';
		Zsub = Zs(idx_phi(j):idx_phi(j+1),idx_theta(i):idx_theta(i+1))';
		
		% IDs of points
		pt_ID = (0:Msub*Nsub-1) + ID_pt_start;
		pt_ID = reshape(pt_ID', Nsub, Msub);
		pt_ID = pt_ID';
		
		% Topology of quadralaterals connecting points
		top_sub1 = pt_ID(1:Msub-1,1:Nsub-1)';
		top_sub2 = pt_ID(1:Msub-1,2:Nsub)';
		top_sub3 = pt_ID(2:Msub,2:Nsub)';
		top_sub4 = pt_ID(2:Msub,1:Nsub-1)';
		
		% Store grid points and topologies of subpatch
		if black
			X_black(ID_pt_start+1:N_pt_black,:) = [Xsub(:), Ysub(:), Zsub(:)];
			topology_black(ID_quad_start+1:N_quad_black,:) = ...
				[top_sub1(:), top_sub2(:), top_sub3(:), top_sub4(:)];
		else
			X_red(ID_pt_start+1:N_pt_red,:) = [Xsub(:), Ysub(:), Zsub(:)];
			topology_red(ID_quad_start+1:N_quad_red,:) = ...
				[top_sub1(:), top_sub2(:), top_sub3(:), top_sub4(:)];
		end
	end
end