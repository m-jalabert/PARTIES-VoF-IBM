function [output] = write_Reader_Vector_Velocity_xmf(path, Nx, Ny, Nz, vector_data)

if isfile(fullfile(path, 'Reader_vector.xmf'))
    delete (fullfile(path, 'Reader_vector.xmf'))
end

Nx_i = Nx - 1;
Ny_i = Ny - 1;
Nz_i = Nz - 1; 
N_vec = 3;

no_vector_files = size(vector_data,1);

file = fullfile(path, strcat('Reader_vector.xmf'));

file_begin = fopen(file, 'w');

fprintf(file_begin,'<?xml version=\"1.0\" ?>\n'); 
fprintf(file_begin,'<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n');
fprintf(file_begin,'<Xdmf Version=\"2.0\">\n');
fprintf(file_begin,'  <Domain>\n\n');

fprintf(file_begin,'    <Topology TopologyType="3DRectMesh" NumberOfElements="%d %d %d"/>\n', Nz_i, Ny_i, Nx_i);
fprintf(file_begin,'    <Geometry GeometryType="VXVYVZ">\n');
fprintf(file_begin,'      <DataItem ItemType="HyperSlab" Dimensions="%d" Type="HyperSlab">\n', Nx_i);
fprintf(file_begin,'        <DataItem Dimensions="3" Format="XML">\n');
fprintf(file_begin,'          0    1    %d\n', Nx_i);
fprintf(file_begin,'        </DataItem>\n');
fprintf(file_begin,'        <DataItem Format="HDF" Dimensions="%d">\n', Nx_i); % Nx
fprintf(file_begin,'          Vector_0.h5:/grid/xc\n');
fprintf(file_begin,'        </DataItem>\n');
fprintf(file_begin,'      </DataItem>\n\n');

fprintf(file_begin,'      <DataItem ItemType="HyperSlab" Dimensions="%d" Type="HyperSlab">\n', Ny_i);
fprintf(file_begin,'        <DataItem Dimensions="3" Format="XML">\n');
fprintf(file_begin,'          0    1    %d\n', Ny_i);
fprintf(file_begin,'        </DataItem>\n');
fprintf(file_begin,'        <DataItem Format="HDF" Dimensions="%d">\n', Ny_i); % Ny
fprintf(file_begin,'          Vector_0.h5:/grid/yc\n');
fprintf(file_begin,'        </DataItem>\n');
fprintf(file_begin,'      </DataItem>\n\n');

fprintf(file_begin,'      <DataItem ItemType="HyperSlab" Dimensions="%d" Type="HyperSlab">\n', Nz_i);
fprintf(file_begin,'        <DataItem Dimensions="3" Format="XML">\n');
fprintf(file_begin,'          0    1    %d\n', Nz_i);
fprintf(file_begin,'        </DataItem>\n');
fprintf(file_begin,'        <DataItem Format="HDF" Dimensions="%d">\n', Nz_i); % Nz
fprintf(file_begin,'          Vector_0.h5:/grid/zc\n');
fprintf(file_begin,'        </DataItem>\n');
fprintf(file_begin,'      </DataItem>\n');
fprintf(file_begin,'    </Geometry>\n\n');

fprintf(file_begin,'    <Grid Name="TemporalGrid" GridType="Collection" CollectionType="Temporal">\n\n');


file_main = fopen(file, 'a');

for i = 1 : no_vector_files

    fprintf(file_main, '      <Grid Name="SpatialGrid_%d" GridType="Uniform">\n', i-1);
    fprintf(file_main, '        <Time Value="%f"/>\n', vector_data(i,2));
    fprintf(file_main, '        <Topology Reference="/Xdmf/Domain/Topology[1]"/>\n');
    fprintf(file_main, '        <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>\n');
    fprintf(file_main, '        <Attribute Name="vector" AttributeType="Vector" Center="Node">\n');
    fprintf(file_main, '          <DataItem ItemType="HyperSlab" Dimensions="%d %d %d %d" Type="HyperSlab">\n', Nz_i, Ny_i, Nx_i, N_vec);
    fprintf(file_main, '            <DataItem Dimensions="3 4" Format="XML">\n');
    fprintf(file_main, '              0    0    0    0 \n');
    fprintf(file_main, '              1    1    1    1 \n');
    fprintf(file_main, '              %-4d %-4d %-4d %-4d\n', Nz_i, Ny_i, Nx_i, N_vec);
    fprintf(file_main, '            </DataItem>\n');
    fprintf(file_main, '            <DataItem Format="HDF" NumberType="Double" Precision="8" Dimensions="%d %d %d %d">\n', Nx_i, Ny_i, Nz_i, N_vec); % Nx, Ny, Nz, N_vec
    fprintf(file_main, '              Vector_%d.h5:/vector_velocity\n', i-1);
    fprintf(file_main, '            </DataItem>\n');
    fprintf(file_main, '          </DataItem>\n');
    fprintf(file_main, '        </Attribute>\n');
    fprintf(file_main, '      </Grid>\n\n');

end % i

file_end = fopen(file, 'a');

fprintf(file_end, '    </Grid>\n');
fprintf(file_end, '  </Domain>\n');
fprintf(file_end, '</Xdmf>\n');



output = 'Wrote  Reader_vector.xmf\n';



end

