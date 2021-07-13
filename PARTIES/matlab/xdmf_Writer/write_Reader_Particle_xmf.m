function [output] = write_Reader_Particle_xmf(path, time, type, Type, Np, att_list)

if isfile(fullfile(path, strcat('Reader_p_', type, '.xmf')))
    delete (fullfile(path, strcat('Reader_p_', type, '.xmf')))
end

no_particle_files = numel(dir(fullfile(path, 'Particle_*')));

file = fullfile(path, strcat('Reader_p_', type, '.xmf'));

file_begin = fopen(file, 'w');

fprintf(file_begin,'<?xml version="1.0" ?>\n');
fprintf(file_begin,'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n');
fprintf(file_begin,'<Xdmf Version="2.0">\n');
fprintf(file_begin,'  <Domain>\n\n');

fprintf(file_begin,'    <Grid Name="TemporalGrid" GridType="Collection" CollectionType="Temporal">\n');

file_main = fopen(file, 'a');

for i = 1 : no_particle_files

    fprintf(file_main, '      <Grid Name="%sGrid_%d">\n', Type, i-1);
    fprintf(file_main, '        <Time Value="%f"/>\n', time(i));
    fprintf(file_main, '        <Topology Type="Polyvertex" NumberOfElements="%d" />\n\n', Np);
    
    fprintf(file_main, '        <Geometry Type="XYZ">\n');
    fprintf(file_main, '          <DataItem Format="HDF" Dimensions="%d 3">\n', Np);
    fprintf(file_main, '            Particle_%d.h5:/%s/X\n', i-1, type);
    fprintf(file_main, '          </DataItem>\n');
    fprintf(file_main, '        </Geometry>\n');
    
    for j = 1 : size(att_list, 1)
        
        fprintf(file_main, '\n        <Attribute Name="%s" AttributeType="%s" Center="Node">\n', att_list{j,1}, att_list{j,3});
        fprintf(file_main, '          <DataItem Format="HDF" NumberType="Double" Dimensions="%d %d">\n', Np, att_list{j,2});
        fprintf(file_main, '            Particle_%d.h5:/%s/%s\n', i-1, type, att_list{j,1});
        fprintf(file_main, '          </DataItem>\n');
        fprintf(file_main, '        </Attribute>\n');
        
        
    end % j
    
    fprintf(file_main, '      </Grid>\n\n');

end % i

file_end = fopen(file, 'a');

fprintf(file_end, '    </Grid>\n');
fprintf(file_end, '  </Domain>\n');
fprintf(file_end, '</Xdmf>\n');

output = strcat('Wrote  Reader_p_', type ,'.xmf\n');

end

