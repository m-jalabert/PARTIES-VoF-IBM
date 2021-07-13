clear;

%% Input
% enter path to simulation folder
% if the whole file "xdmf_Writer" has been copied to the simulation folder,
% set: path = '..';
path = '..'; 


%% Pre-Check
% Check if entered path is corrrect by checking whether Data_0.h5 exists 
if isfile(fullfile(path, 'Data_0.h5'))
    %% Store and print Nx, Ny, and Nz
    Nx = h5read(fullfile(path, 'Data_0.h5'), '/grid/NX');
    Ny = h5read(fullfile(path, 'Data_0.h5'), '/grid/NY');
    Nz = h5read(fullfile(path, 'Data_0.h5'), '/grid/NZ');

    fprintf('##########################\n')
    fprintf('NX = %d \n', Nx)
    fprintf('NY = %d \n', Ny)
    fprintf('NZ = %d \n', Nz)
    fprintf('\n')


    %% Store and print time steps
    no_data_files = numel(dir(fullfile(path, 'Data_*')));
    time = zeros(no_data_files,1);

    for i = 1 : no_data_files
        data_file = strcat('Data_', string(i-1), '.h5');
        % store time steps in an array to provide all time steps to the
        % functions "write_Reader_*"
        time(i) = double(h5read(fullfile(path, data_file), '/time'));
        fprintf('t_%d = %.6f \n', i-1, time(i))
    end %i
    fprintf('\n')


    %% Run "write_Reader_*" functions and print done jobs
    particle_file = h5info(fullfile(path, 'Particle_0.h5')); 
    no_datasets = size(particle_file.Groups,1);

    if no_datasets > 1
    
        % loop through all particle-datasets
        % important when at least one fixed particle and one mobile particle exist
        for j = 2 : no_datasets
            % Check if particle is fixed or mobile
            if particle_file.Groups(j).Name == "/fixed"
                type = 'fixed';
                Type = 'Fixed';
            else
                type = 'mobile';
                Type = 'Mobile';
            end

            % Number of particles
            Np = size(h5read(fullfile(path, 'Particle_0.h5'), strcat('/', type, '/R')), 2);

            fprintf('Np_%s = %d \n', type, Np)

            size_att = size(particle_file.Groups(j).Datasets, 1) - 1;
            att_list = cell(size_att, 3);

            % store data in an attribute list with which all information
            % are provided to the function "write_Reader_Particle.xmf"
            for l = 1 : size_att
                att_list{l,1} = particle_file.Groups(j).Datasets(l).Name;
                att_list{l,2} = particle_file.Groups(j).Datasets(l).Dataspace.Size(1);

                if att_list{l,2} == 1
                    att_list{l,3} = 'Scalar';
                elseif att_list{l,2} == 3
                    att_list{l,3} = 'Vector';
                elseif att_list{l,2} == 9
                    att_list{l,3} = 'Tensor';
                end

            end %k

            % write Reader_p_fixed/mobile.xmf
            p = write_Reader_Particle_xmf(path, time, type, Type, Np, att_list);
            fprintf(p)
            fprintf('\n')

        end %j

    end % if
    
    
    % write Reader_c.xmf
    c = write_Reader_c_xmf(path, Nx, Ny, Nz, time);
    fprintf(c)

    % write Reader_u.xmf
    u = write_Reader_Scalar_Velocity_xmf(path, Nx, Ny, Nz, time, "u");
    fprintf(u)

    % write Reader_v.xmf
    v = write_Reader_Scalar_Velocity_xmf(path, Nx, Ny, Nz, time, "v");
    fprintf(v)

    % write Reader_w.xmf
    w = write_Reader_Scalar_Velocity_xmf(path, Nx, Ny, Nz, time, "w");
    fprintf(w)

    % check if Vector_*.h5 files exist
    vector_dir = dir(fullfile(path, 'Vector_*'));
    size_vector = size(vector_dir,1);
        
    if size_vector > 0
        % write Reader_vector.xmf
        vector = write_Reader_Vector_Velocity_xmf(path, Nx, Ny, Nz, time);
        fprintf(vector)
    end

    fprintf('##########################\n')

else
    
    fprintf('##########################\n')
    fprintf('Please check the entered path, because no simulation data (h5 files) could be found!\n')
    fprintf('##########################\n')

end
