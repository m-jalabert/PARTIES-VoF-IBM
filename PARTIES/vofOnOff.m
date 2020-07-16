clear all;
close all;

fld = {'VOFon'};

ns = 0;
ne  = 14;
for i  = 1:length(fld)
    
    wd = fld{i};
    cnt = 0;
    
    for n = ns:ne
    cnt = cnt + 1;
    c0 = h5read([wd '/Data_' num2str(n) '.h5'],'/Conc/0');
    x = h5read([wd '/Data_' num2str(n) '.h5'],'/grid/xc');
    y = h5read([wd '/Data_' num2str(n) '.h5'],'/grid/yc');
    vf = h5read([wd '/Data_' num2str(n) '.h5'],'/vfc');
    mass(cnt) = sum(sum(sum(c0.*(1-vf))));
    t(cnt) =  h5read([wd '/Data_' num2str(n) '.h5'],'/time');
    end
    
    figure(1)
    plot(t,(mass-mass(1))/mass(1)); hold on;
    
    figure
    h = pcolor(x,y,c0(:,:,10)');
    figure
    h = pcolor(x,y,vf(:,:,10)');
    
end