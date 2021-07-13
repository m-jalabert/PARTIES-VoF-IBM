function [output] = write_Reader_c_xmf(path, Nx, Ny, Nz, time)


if isfile(fullfile(path, 'Reader_c.xmf'))
    delete (fullfile(path, 'Reader_c.xmf'))
end

Nx_i = Nx - 1;
Ny_i = Ny - 1;
Nz_i = Nz - 1;   


no_data_files = numel(dir(fullfile(path, 'Data_*')));

file = fullfile(path, 'Reader_c.xmf');

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
fprintf(file_begin,'          Data_0.h5:/grid/xc\n');
fprintf(file_begin,'        </DataItem>\n');
fprintf(file_begin,'      </DataItem>\n\n');

fprintf(file_begin,'      <DataItem ItemType="HyperSlab" Dimensions="%d" Type="HyperSlab">\n', Ny_i);
fprintf(file_begin,'        <DataItem Dimensions="3" Format="XML">\n');
fprintf(file_begin,'          0    1    %d\n', Ny_i);
fprintf(file_begin,'        </DataItem>\n');
fprintf(file_begin,'        <DataItem Format="HDF" Dimensions="%d">\n', Ny);
fprintf(file_begin,'          Data_0.h5:/grid/yc\n');
fprintf(file_begin,'        </DataItem>\n');
fprintf(file_begin,'      </DataItem>\n\n');

fprintf(file_begin,'      <DataItem ItemType="HyperSlab" Dimensions="%d" Type="HyperSlab">\n', Nz_i);
fprintf(file_begin,'        <DataItem Dimensions="3" Format="XML">\n');
fprintf(file_begin,'          0    1    %d\n', Nz_i);
fprintf(file_begin,'        </DataItem>\n');
fprintf(file_begin,'        <DataItem Format="HDF" Dimensions="%d">\n', Nz);
fprintf(file_begin,'          Data_0.h5:/grid/zc\n');
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
    fprintf(file_main, '        <Attribute Name="p" AttributeType="Scalar" Center="Node">\n');
    fprintf(file_main, '          <DataItem ItemType="HyperSlab" Dimensions="%d %d %d" Type="HyperSlab">\n', Nz_i, Ny_i, Nx_i);
    fprintf(file_main, '            <DataItem Dimensions="3 3" Format="XML">\n');
    fprintf(file_main, '              0    0    0 \n');
    fprintf(file_main, '              1    1    1\n');
    fprintf(file_main, '              %-4d %-4d %-4d\n', Nz_i, Ny_i, Nx_i);
    fprintf(file_main, '            </DataItem>\n');
    fprintf(file_main, '            <DataItem Format="HDF" NumberType="Double" Precision="8" Dimensions="%d %d %d">\n', Nx, Ny, Nz);
    fprintf(file_main, '              Data_%d.h5:/p\n', i-1);
    fprintf(file_main, '            </DataItem>\n');
    fprintf(file_main, '          </DataItem>\n');
    fprintf(file_main, '        </Attribute>\n');
    fprintf(file_main, '      </Grid>\n\n');

end % i

file_end = fopen(file, 'a');

fprintf(file_end, '    </Grid>\n');
fprintf(file_end, '  </Domain>\n');
fprintf(file_end, '</Xdmf>\n');



output = 'Wrote  Reader_c.xmf\n';

end

