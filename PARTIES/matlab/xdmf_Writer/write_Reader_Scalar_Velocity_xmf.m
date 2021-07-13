function [output] = write_Reader_Scalar_Velocity_xmf(path, Nx, Ny, Nz, time, var)


if isfile(fullfile(path, strcat('Reader_', var, '.xmf')))
    delete (fullfile(path, strcat('Reader_', var, '.xmf')))
end

if var == "u"
    N = Nx;
    Nx_i = Nx;
    Ny_i = Ny - 1;
    Nz_i = Nz - 1;   
    x = 'xu';
    y = 'yc';
    z = 'zc';
elseif var == "v"
    N = Ny;
    Nx_i = Nx - 1;
    Ny_i = Ny;
    Nz_i = Nz - 1; 
    x = 'xc';
    y = 'yv';
    z = 'zc';
elseif var == "w"
    N = Nz;
    Nx_i = Nx - 1;
    Ny_i = Ny - 1;
    Nz_i = Nz; 
    x = 'xc';
    y = 'yc';
    z = 'zw';
end


no_data_files = numel(dir(fullfile(path, 'Data_*')));

file = fullfile(path, strcat('Reader_', var, '.xmf'));

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
fprintf(file_begin,'        <DataItem Format="HDF" Dimensions="%d">\n', Nx);
fprintf(file_begin,'          Data_0.h5:/grid/%s\n', x);
fprintf(file_begin,'        </DataItem>\n');
fprintf(file_begin,'      </DataItem>\n\n');

fprintf(file_begin,'      <DataItem ItemType="HyperSlab" Dimensions="%d" Type="HyperSlab">\n', Ny_i);
fprintf(file_begin,'        <DataItem Dimensions="3" Format="XML">\n');
fprintf(file_begin,'          0    1    %d\n', Ny_i);
fprintf(file_begin,'        </DataItem>\n');
fprintf(file_begin,'        <DataItem Format="HDF" Dimensions="%d">\n', Ny);
fprintf(file_begin,'          Data_0.h5:/grid/%s\n', y);
fprintf(file_begin,'        </DataItem>\n');
fprintf(file_begin,'      </DataItem>\n\n');

fprintf(file_begin,'      <DataItem ItemType="HyperSlab" Dimensions="%d" Type="HyperSlab">\n', Nz_i);
fprintf(file_begin,'        <DataItem Dimensions="3" Format="XML">\n');
fprintf(file_begin,'          0    1    %d\n', Nz_i);
fprintf(file_begin,'        </DataItem>\n');
fprintf(file_begin,'        <DataItem Format="HDF" Dimensions="%d">\n', Nz);
fprintf(file_begin,'          Data_0.h5:/grid/%s\n', z);
fprintf(file_begin,'        </DataItem>\n');
fprintf(file_begin,'      </DataItem>\n');
fprintf(file_begin,'    </Geometry>\n\n');

fprintf(file_begin,'    <Grid Name="TemporalGrid" GridType="Collection" CollectionType="Temporal">\n\n');


file_main = fopen(file, 'a');

for i = 1 : no_data_files

    fprintf(file_main, '      <Grid Name="SpatialGrid_%d" GridType="Uniform">\n', i-1);
    fprintf(file_main, '        <Time Value="%f"/>\n', time(i));
    fprintf(file_main, '        <Topology Reference="/Xdmf/Domain/Topology[1]"/>\n');
    fprintf(file_main, '        <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>\n');
    fprintf(file_main, '        <Attribute Name="%s" AttributeType="Scalar" Center="Node">\n', var);
    fprintf(file_main, '          <DataItem ItemType="HyperSlab" Dimensions="%d %d %d" Type="HyperSlab">\n', Nz_i, Ny_i, Nx_i);
    fprintf(file_main, '            <DataItem Dimensions="3 3" Format="XML">\n');
    fprintf(file_main, '              0    0    0 \n');
    fprintf(file_main, '              1    1    1\n');
    fprintf(file_main, '              %-4d %-4d %-4d\n', Nz_i, Ny_i, Nx_i);
    fprintf(file_main, '            </DataItem>\n');
    fprintf(file_main, '            <DataItem Format="HDF" NumberType="Double" Precision="8" Dimensions="%d %d %d">\n', Nx, Ny, Nz);
    fprintf(file_main, '              Data_%d.h5:/%s\n', i-1, var);
    fprintf(file_main, '            </DataItem>\n');
    fprintf(file_main, '          </DataItem>\n');
    fprintf(file_main, '        </Attribute>\n');
    fprintf(file_main, '        <Attribute Name="vf%s" AttributeType="Scalar" Center="Node">\n', var);
    fprintf(file_main, '          <DataItem ItemType="HyperSlab" Dimensions="%d %d %d" Type="HyperSlab">\n', Nz_i, Ny_i, Nx_i);
    fprintf(file_main, '            <DataItem Dimensions="3 3" Format="XML">\n');
    fprintf(file_main, '              0    0    0 \n');
    fprintf(file_main, '              1    1    1\n');
    fprintf(file_main, '              %-4d %-4d %-4d\n', Nz_i, Ny_i, Nx_i);
    fprintf(file_main, '            </DataItem>\n');
    fprintf(file_main, '            <DataItem Format="HDF" NumberType="Double" Precision="8" Dimensions="%d %d %d">\n', Nx, Ny, Nz);
    fprintf(file_main, '              Data_%d.h5:/vf%s\n', i-1, var);
    fprintf(file_main, '            </DataItem>\n');
    fprintf(file_main, '          </DataItem>\n');
    fprintf(file_main, '        </Attribute>\n');
    fprintf(file_main, '      </Grid>\n\n');

end % i

file_end = fopen(file, 'a');

fprintf(file_end, '    </Grid>\n');
fprintf(file_end, '  </Domain>\n');
fprintf(file_end, '</Xdmf>\n');





output = strcat('Wrote  Reader_', var ,'.xmf\n');

end

