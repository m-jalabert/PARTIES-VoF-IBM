function [] = add_periodic_particles(hdf5_idx, Nproc)
% ADD_PERIODIC_PARTICLES   Duplicates periodic particles in HDF5 file
%    ADD_PERIODIC_PARTICLES(hdf5_idx, Nproc) Edits all Particle_*.h5 files
%    where * is a vector of integers specified by the input 'hdf5_idx'.
%    'Nproc' is the number of processors to run on, to speed up the work.
%
%    A new group, 'ID', is also added to the HDF5 file, which contains
%    information about the particle IDs, including for the newly-added particles
%
%    NOTE: y-periodicity is currently unsupported
%
%    Example:
%       add_periodic_particles(0:50, 2)
%
%    duplicates periodic particles for HDF5 files Particle_0.h5 through
%    Particle_50.h5 using two processors.

if (Nproc > 1)
	poolobj = parpool(Nproc);
	parfor count = 1:length(hdf5_idx)
		add_periodic_particles_serial(hdf5_idx(count))
	end
	delete(poolobj)
else
	for count = 1:length(hdf5_idx)
		add_periodic_particles_serial(hdf5_idx(count))
	end
end


% ******************************************************************************
% Add periodic particles to file indicated by 'hdf5_idx'
% ******************************************************************************
function [] = add_periodic_particles_serial(hdf5_idx)

filename = sprintf('%s/Particle_%d.h5', pwd, hdf5_idx);
info = h5info(filename,'/');
time = h5read(filename, '/time');

% Read in domain information
domain.xmin = h5read(filename, '/domain/xmin');
domain.xmax = h5read(filename, '/domain/xmax');
domain.ymin = h5read(filename, '/domain/ymin');
domain.ymax = h5read(filename, '/domain/ymax');
domain.zmin = h5read(filename, '/domain/zmin');
domain.zmax = h5read(filename, '/domain/zmax');
domain.periodic = h5read(filename, '/domain/periodic');

% Find number of fixed and mobile particles
Npf = 0;
Npm = 0;
for j=1:length(info.Groups)
	group = info.Groups(j).Name;
	
	if (strcmp(group, '/mobile') || strcmp(group, '/fixed'))
		
		% Check if file has been edited already
		ginfo = h5info(filename, group);
		for i=1:length(ginfo.Datasets)
			if strcmp(ginfo.Datasets(i).Name, 'ID')
				fprintf('%s already has periodic particles! Skipping...\n', filename)
				return
			end
		end
		
		dsinfo = h5info(filename, [group '/R']);
		if strcmp(group, '/mobile')
			Npm = dsinfo.Dataspace.Size(2);
		elseif strcmp(group, '/fixed')
			Npf = dsinfo.Dataspace.Size(2);
		end
	end
end

tempfile = sprintf('%s/temp_%d.h5', pwd, hdf5_idx);
h5create(tempfile,'/time',1);
h5write(tempfile,'/time',time)

fprintf('Editing %s\n', filename)

for j=1:length(info.Groups)

	group = info.Groups(j).Name;
	if (strcmp(group, '/mobile') || strcmp(group, '/fixed'))

		R = h5read(filename, [group '/R'])';
		X = h5read(filename, [group '/X'])';
		ID = (1:length(R))';

		R_new = [];
		X_new = [];
		ID_new = [];

		if (domain.periodic(1))
			% Particles that extend past lower x domain
			X_temp_idx = find(X(:,1) - R < domain.xmin);
			X_temp = X(X_temp_idx,:);
			R_temp = R(X_temp_idx);
			X_temp_ID = ID(X_temp_idx);

			X_temp(:,1) = X_temp(:,1) + (domain.xmax - domain.xmin);
			X_new = [X_new; X_temp];
			R_new = [R_new; R_temp];
			ID_new = [ID_new; X_temp_ID];

			% Particles that extend past upper x domain
			X_temp_idx = find(X(:,1) + R > domain.xmax);
			X_temp = X(X_temp_idx,:);
			R_temp = R(X_temp_idx);
			X_temp_ID = ID(X_temp_idx);

			X_temp(:,1) = X_temp(:,1) - (domain.xmax - domain.xmin);
			X_new = [X_new; X_temp];
			R_new = [R_new; R_temp];
			ID_new = [ID_new; X_temp_ID];
		end

		% Have x-periodic particles operated on for z-periodic as well
		X = [X; X_new];
		R = [R; R_new];
		ID = [ID; ID_new];

		if (domain.periodic(3))
			% Particles that extend past lower z domain
			X_temp_idx = find(X(:,3) - R < domain.zmin);
			X_temp = X(X_temp_idx,:);
			X_temp_ID = ID(X_temp_idx);

			X_temp(:,3) = X_temp(:,3) + (domain.zmax - domain.zmin);
			X_new = [X_new; X_temp];
			ID_new = [ID_new; X_temp_ID];

			% Particles that extend past upper z domain
			X_temp_idx = find(X(:,3) + R > domain.zmax);
			X_temp = X(X_temp_idx,:);
			X_temp_ID = ID(X_temp_idx);

			X_temp(:,3) = X_temp(:,3) - (domain.zmax - domain.zmin);
			X_new = [X_new; X_temp];
			ID_new = [ID_new; X_temp_ID];
		end

		add_new_particles(filename, tempfile, group, X_new, ID_new);
		add_particle_IDs(tempfile, group, Npf, Npm, ID_new);
	end
end

delete(filename);
movefile(tempfile,filename);


% ******************************************************************************
% Appends periodic particles to master particles in HDF5 file
% ******************************************************************************
function [] = add_new_particles(filename, tempfile, group, X_new, X_new_ID)

info = h5info(filename, group);
for i=1:length(info.Datasets)
	
	datasetname = info.Datasets(i).Name;
	dataset = [group '/' datasetname];
	
	data = h5read(filename, dataset)';
	
	if (strcmp(datasetname, 'X'))
		data = [data; X_new];
	else
		data = [data; data(X_new_ID,:)];
	end
	data = data';
	
	h5create(tempfile, dataset, size(data));
	h5write(tempfile, dataset, data);
	
end


% ******************************************************************************
% Add data structure 'ID' to HDF5 file
% ******************************************************************************
function [] = add_particle_IDs(tempfile, group, Npf, Npm, ID_new)

if strcmp(group, '/fixed')
	ID = (0:Npf-1)';
elseif strcmp(group, '/mobile')
	ID = (Npf:(Npm+Npf-1))';
end

ID = [ID; ID(ID_new)];
ID = ID';

dataset = [group '/ID'];
h5create(tempfile, dataset, size(ID), 'Datatype', 'int32');
h5write(tempfile, dataset, int32(ID));

