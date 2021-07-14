clear;

%% Input 
% enter path to simulation folder
% if matlab script has been copied to simulation folder, set: path = '';
path = '/home/Desktop/PARTIES';


%% Create Vector_*.h5 files
data_dir = dir(fullfile(path, 'Data_*'));
no_files = size(data_dir,1);


if no_files >= 1
    
    fprintf('##########################\n')

    for id = 0 : (no_files-1)
        %% Delete Vector-files if existing
        if isfile(fullfile(path, strcat('Vector_', string(id), '.h5')))
            delete (fullfile(path, strcat('Vector_', string(id), '.h5')))
        end

        %% Read data of existing h5 files
        data_file = strcat('Data_', string(id), '.h5');

        u = h5read(fullfile(path, data_file), '/u');
        v = h5read(fullfile(path, data_file), '/v');
        w = h5read(fullfile(path, data_file), '/w');

        xc = h5read(fullfile(path, data_file), '/grid/xc');
        yc = h5read(fullfile(path, data_file), '/grid/yc');
        zc = h5read(fullfile(path, data_file), '/grid/zc');

        Nx = h5read(fullfile(path, data_file), '/grid/NX');
        Ny = h5read(fullfile(path, data_file), '/grid/NY');
        Nz = h5read(fullfile(path, data_file), '/grid/NZ');


        %% Create hdf5 file
        file = strcat('Vector_', string(id), '.h5');
        create_file = H5F.create(fullfile(path, file));


        %% Dataset time
        time = double(h5read(fullfile(path, data_file), '/time'));
        dim_time = [1];
        space_time = H5S.create_simple(1, dim_time, dim_time);
        dataset_time = H5D.create(create_file,'/time','H5T_NATIVE_DOUBLE',space_time,'H5P_DEFAULT');

        H5D.write(dataset_time, 'H5ML_DEFAULT','H5S_ALL','H5S_ALL', 'H5P_DEFAULT', time);


        %% Dataset velocity vector
        % The velocity vector is cell-centered
        Nx_i = Nx - 1;
        Ny_i = Ny - 1;
        Nz_i = Nz - 1;

        dims_vector = [double(Nz_i) double(Ny_i) double(Nx_i) 3]; % Z, Y, X, vector
        space_vector = H5S.create_simple(4,dims_vector,dims_vector);
        dataset_vector = H5D.create(create_file,'/vector_velocity','H5T_NATIVE_DOUBLE',space_vector,'H5P_DEFAULT');

        velocity_matrix_4D = zeros (3, Nx_i, Ny_i, Nz_i); % create empty 4D matrix to store velocity components
        % structure is vice versa to h5-file

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

        H5D.write(dataset_vector, 'H5ML_DEFAULT','H5S_ALL','H5S_ALL', 'H5P_DEFAULT', velocity_matrix_4D);


        %% Create group "grid"
        group = H5G.create(create_file,'grid','H5P_DEFAULT','H5P_DEFAULT','H5P_DEFAULT');


        %% Datasets /grid/Nx,Ny,Nz
        dim_N = [1];

        space_Nx = H5S.create_simple(1, dim_N, dim_N);
        dataset_Nx = H5D.create(create_file,'/grid/NX','H5T_NATIVE_INT',space_Nx,'H5P_DEFAULT');
        H5D.write(dataset_Nx, 'H5ML_DEFAULT','H5S_ALL','H5S_ALL', 'H5P_DEFAULT', Nx);

        space_Ny = H5S.create_simple(1, dim_N, dim_N);
        dataset_Ny = H5D.create(create_file,'/grid/NY','H5T_NATIVE_INT',space_Ny,'H5P_DEFAULT');
        H5D.write(dataset_Ny, 'H5ML_DEFAULT','H5S_ALL','H5S_ALL', 'H5P_DEFAULT', Ny);

        space_Nz = H5S.create_simple(1, dim_N, dim_N);
        dataset_Nz= H5D.create(create_file,'/grid/NZ','H5T_NATIVE_INT',space_Nz,'H5P_DEFAULT');
        H5D.write(dataset_Nz, 'H5ML_DEFAULT','H5S_ALL','H5S_ALL', 'H5P_DEFAULT', Nz);


        %% Datasets /grid/xc,yc,zc
        xci = xc;
        xci(end) = [];
        dim_xci = [size(xci,1)];
        space_xci = H5S.create_simple(1, dim_xci, dim_xci);
        dataset_xci = H5D.create(create_file,'/grid/xc','H5T_NATIVE_DOUBLE',space_xci,'H5P_DEFAULT');
        H5D.write(dataset_xci, 'H5ML_DEFAULT','H5S_ALL','H5S_ALL', 'H5P_DEFAULT', double(xci));

        yci = yc;
        yci(end) = [];
        dim_yci = [size(yci,1)];
        space_yci = H5S.create_simple(1, dim_yci, dim_yci);
        dataset_yci = H5D.create(create_file,'/grid/yc','H5T_NATIVE_DOUBLE',space_yci,'H5P_DEFAULT');
        H5D.write(dataset_yci, 'H5ML_DEFAULT','H5S_ALL','H5S_ALL', 'H5P_DEFAULT', double(yci));

        zci = zc;
        zci(end) = [];
        dim_zci = [size(zci,1)];
        space_zci = H5S.create_simple(1, dim_zci, dim_zci);
        dataset_zci = H5D.create(create_file,'/grid/zc','H5T_NATIVE_DOUBLE',space_zci,'H5P_DEFAULT');
        H5D.write(dataset_zci, 'H5ML_DEFAULT','H5S_ALL','H5S_ALL', 'H5P_DEFAULT', double(zci));


        %% Close h5 libraries
        H5D.close(dataset_Nx);
        H5D.close(dataset_Ny);
        H5D.close(dataset_Nz);
        H5D.close(dataset_time);
        H5D.close(dataset_vector);
        H5D.close(dataset_xci);
        H5D.close(dataset_yci);
        H5D.close(dataset_zci);

        H5F.close(create_file);

        H5G.close(group);

        H5S.close(space_Nx);
        H5S.close(space_Ny);
        H5S.close(space_Nz);
        H5S.close(space_time);
        H5S.close(space_vector);
        H5S.close(space_xci);
        H5S.close(space_yci);
        H5S.close(space_zci);


        fprintf('%s done\n', file);

    end %id

    fprintf('##########################\n')

else
    fprintf('##########################\n')
    fprintf('Please check the entered path, because no Data_*.h5 files could be found!\n')
    fprintf('##########################\n')
end
    


