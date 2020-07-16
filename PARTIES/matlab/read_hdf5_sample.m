%%% Diffusion with a sphere in the center fro VOF development

clear;
addpath '/home/thko6956/2017/PARTIES/Test-Particle/matlab_scrips/'


path.case1= '/home/thko6956/2017/PARTIES/Test-Particle/partest/';   %% to be chnaged
path.case2= '/home/thko6956/2017/PARTIES/Test-Particle/cases/Test_swimmers_stratified/original/'




DO3D = 1;
DO2D = 0;

distribution_2d=0



caases={ 'case1', 'case2'};  


for idc=1:length(caases) 
    
caase= caases{idc}
datapath=path.(caase);



if (DO3D == 1)
    
filelist= dir([datapath 'Data_*.h5'] );

if (~isempty(filelist))

n_timesteps=length(filelist)
num=1 %%% read only every num th file    
       
integral_conc_fluid.(caase)=[];
integral_conc_full.(caase)=[];
time.(caase)=[];
kine_f.(caase)=[];
kine.(caase)=[]
dyc_bot.(caase)=[];
mass_p.(caase)=[];
mass_f.(caase)=[];
mass.(caase)=[];
currentfile=[datapath filelist(1).name];
xu.(caase)=   h5read(currentfile,['/grid/xu']); 
xc.(caase)=   h5read(currentfile,['/grid/xc']);   
yv.(caase)=   h5read(currentfile,['/grid/yv']); 
yc.(caase)=   h5read(currentfile,['/grid/yc']); 
zw.(caase)=   h5read(currentfile,['/grid/zw']); 
zc.(caase)=   h5read(currentfile,['/grid/zc']); 

Lx=xu.(caase)(end);
Ly=yv.(caase)(end);
Lz=zw.(caase)(end);
Ny.(caase)=length(xu.(caase))-1;


    for index_time=1:length(filelist)/num
        
        currentfile=[datapath filelist(index_time*num).name];
       % currentinfo=h5info(currentfile);
        i=index_time;
         time.(caase)(i)= h5read(currentfile,['/time']);
    end
    
    [ time.(caase) timeindex.(caase)]=sort( time.(caase));
         
    
    for index_time=1:length(filelist)/num
        
          i=timeindex.(caase)(index_time)   
             
        currentfile=[datapath filelist(i*num).name];
       p.(caase){i} = h5read(currentfile,['/p']);
   %     u.(caase){i} = h5read(currentfile,['/u']);
   %     v.(caase){i} = h5read(currentfile,['/v']);
   %     w.(caase){i} = h5read(currentfile,['/w']);
     %   c.(caase){i} = h5read(currentfile,['/Conc/0']);
 
       
    
        
  


     end
    end    
end



    
 
if (DO2D == 1)

filelist= dir([datapath 'Data2d_*.h5'] );

if (length(filelist) > 0)

n_timesteps=length(filelist)
num=1 %%% read only every num th file    
       


currentfile=[datapath filelist(1).name];
currentinfo=h5info(currentfile);

xc.(caase)=   h5read(currentfile,['/grid/x']); 
yc.(caase)=   h5read(currentfile,['/grid/y']); 
zc.(caase)=   h5read(currentfile,['/grid/z']); 

Lx=0.5*(xc.(caase)(end)+xc.(caase)(end-1));
Ly=0.5*(yc.(caase)(end)+yc.(caase)(end-1));
Lz=0.5*(zc.(caase)(end)+yc.(caase)(end-1));
Nx.(caase)=length(xc.(caase))-1;
Ny.(caase)=length(yc.(caase))-1;
Nz.(caase)=length(zc.(caase))-1;
time2d.(caase) = [];
frontlocation.(caase)=[];
frontlocation2.(caase)=[];
integral_kinEnergy_fluid.(caase)=[];
integral_kinEnergy_part.(caase)=[];
integral_nusselt_c0.(caase)=[];
integral_nusselt_c1.(caase)=[];
integral_potEnergy_fluid_c0.(caase)=[];
integral_potEnergy_fluid_c1.(caase)=[];
integral_potEnergy_part.(caase)=[];
integral_viscdiss.(caase)= [];
integral_buoyantwork.(caase)= [];
viscdiss.(caase)=[];
integral_buoyantwork_c0.(caase)=[];
integral_buoyantwork_c1.(caase)=[];
integral_grad_y_c0.(caase)=[];
grad_y_c0.(caase)=[];
integral_grad_y_c0.(caase)=[];
integral_flux_uyc0.(caase)=[];
integral_flux_uyc1.(caase)=[];
flux_uyc0.(caase)=[];
flux_uyc1.(caase)=[];
front_speed.(caase)=[];
    for i=1:length(filelist)/num
      currentfile=[datapath filelist(i*num).name];  
      time2d.(caase)(i)= h5read(currentfile,['/time']);
    end
    
    [tmp , readindex]=sort(time2d.(caase));
    for index_time=1:length(filelist)/num
        
        ii=readindex(index_time);
        currentfile=[datapath filelist(ii*num).name];
        i=index_time;
        current_height_full.(caase){i} = h5read(currentfile,['/xdim/Current_height_fluid_c0']);
        current_height_fluid.(caase){i} = h5read(currentfile,['/xdim/Current_height_full_c0']);
        
        if(distribution_2d == 1)
        conc0_mean.(caase)(:,:,i) = h5read(currentfile,['/c0Mean']);
         viscDiss.(caase)(:,:,i) = h5read(currentfile,['/viscDiss']);
       % grad_y_c0.(caase)(:,:,i) = h5read(currentfile,['/grad_y_c0']);
        flux_uyc0.(caase)(:,:,i) = h5read(currentfile,['/flux_vc0']);
  %      flux_uyc1.(caase)(:,:,i) = h5read(currentfile,['/flux_vc1']);
   %     fluid_vf.(caase)(:,:,i) = h5read(currentfile,['/vf_avg']);
       %  conc0_clice.(caase)(:,:,i) = h5read(currentfile,['/c0Slice']);
     %   integral_grad_y_c0.(caase)(i) = mean(mean(grad_y_c0.(caase)(1:end-1,1:end-1,i)));
      %  integral_flux_uyc0.(caase)(i)=mean(mean(flux_uyc0.(caase)(1:end-1,1:end-1,i).*fluid_vf.(caase)(1:end-1,1:end-1,i) ));
      %  integral_flux_uyc1.(caase)(i)=mean(mean(flux_uyc1.(caase)(1:end-1,1:end-1,i).*fluid_vf.(caase)(1:end-1,1:end-1,i) ))
        end
        
        integral_kinEnergy_fluid.(caase)(i)=h5read(currentfile,['/integral/integral_kinEnergy_fluid']);
        integral_nusselt_c0.(caase)(i)=h5read(currentfile,['/integral/integral_nusselt_c0']);
        %integral_nusselt_c1.(caase)=h5read(currentfile,['/integral/integral_nusselt_c1']);
        integral_potEnergy_fluid_c0.(caase)(i)=h5read(currentfile,['/integral/integral_potEnergy_fluid_c0']);
        %integral_potEnergy_fluid_c1.(caase)(i)=h5read(currentfile,['/integral/integral_potEnergy_fluid_c1']);
        integral_viscdiss.(caase)(i)=h5read(currentfile, ['/integral/integral_viscDiss']);
        integral_buoyantwork_c0.(caase)(i)=h5read(currentfile, ['/integral/integral_buoyantwork_c0']);
    %    integral_buoyantwork_c1.(caase)(i)=h5read(currentfile, ['/integral/integral_buoyantwork_c1']);
        integral_conc_fluid_c0.(caase)(i)=h5read(currentfile,['/integral/integral_conc_fluid_c0']);
     %   integral_conc_fluid_c1.(caase)(i)=h5read(currentfile,['/integral/integral_conc_fluid_c1']);
        integral_conc_full_c0.(caase)(i)=h5read(currentfile,['/integral/integral_conc_full_c0']);
      %  integral_conc_full_c1.(caase)(i)=h5read(currentfile,['/integral/integral_conc_full_c1']);
        
     %   front_location.(caase)(i) = h5read(currentfile,['/integral/Front_location_c0']);
        time2d.(caase)(i) = h5read(currentfile,['/time']);
          % extra cal
        
     %
end

end
end
    
    
    
    
    
   
end


%%





lineStyles = linspecer(6);
hfig=figure;
width = 3;     % Width in inches
height = 2.5;     % Height in inches


set(0,'defaulttextinterpreter', 'Latex')
fig = gcf;    
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1)*1.1 pos(2)*1.1 width*100, height*100]); %<- Set size
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/1.8;
bottom = (papersize(2)- height)/1.8;
myfiguresize = [left, bottom, width*1.1, height*1.1];
set(gcf,'PaperPosition', myfiguresize);


plot( time2d.caase1,  integral_kinEnergy_fluid.caase1  , 'xk','LineWidth', 1, 'MarkerSize', 5,'Color',lineStyles(1,:));
%hold on
%plot( time2d,   data2 ,'^k','LineWidth', 1, 'MarkerSize', 5,'Color',lineStyles(2,:));

xlabel('$t$','Fontsize',8, 'Interpreter', 'Latex');
ylabel('$Data$','Fontsize',8, 'Interpreter', 'Latex');
set(gca,'Fontsize',8,'YGrid', 'on', 'XGrid', 'on', 'XLim', [0  10], 'YLim', [0, 4]);%, 'XLim', [0 0.4], 'YLim', [0 400]);%, 'YScale', 'Log', 'YLim', [1e-12 1e-5]);%, 'YTick', [1e-10 1e-9 1e-8 1e-7 1e-5 1e-4]);
h=legend('data1', 'data2' );
set(h,'Fontsize',8,'Location','NorthWest', 'Box', 'on', 'Interpreter', 'Latex' );
matlabfrag( [picturepath 'xx'], 'epspad', [10 ,10,0,0]);
%%print(  [picturepath 'impermeable_sphere/diff_comp'],'-depsc','-painters', '-loose' )

%close(fig);



