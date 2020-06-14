%%%Short testing

%%%%%%% Time dependent diffusion test case
%%%%% Processing test cases with matlab
clear

addpath '/home/t.koellner/2017/Test-Particle/matlab_scrips/'
path.heat_exp= '/home/t.koellner/2017/Test-Particle/partest/';   %% to be chnaged
picturepath ='/home/t.koellner/2017/Test-Particle/picture/';




%% First shear flow
%caases={'heat_exp'}
caases={'heat_exp'}
Pe=1
caase= caases{1}
CONC=1;

%pd = h5read([path.(caase) 'pressure'], '/data');
%Vp = h5read([path.(caase) 'velopro'], '/data');
%conc = h5read([path.(caase) 'conc'], '/data');
%Vunp = h5read([path.(caase) 'velounpro'], '/data');



for idc=1:length(caases) 
    
caase= caases{idc}

datapath=path.(caase);
what2read={'p', 'u', 'v','w', 'c'}
% read in all data


filelist= dir([datapath 'Data*.h5'] );
n_timesteps=length(filelist)
num=1  %%% read only every num th file    
        
    for index_time=1:length(filelist)/num
        
        currentfile=[datapath filelist(index_time*num).name];
        currentinfo=h5info(currentfile);
        i=index_time;
        p.(caase){i} = h5read(currentfile,['/p']);
        u.(caase){i} = h5read(currentfile,['/u']);
        v.(caase){i} = h5read(currentfile,['/v']);
        w.(caase){i} = h5read(currentfile,['/w']);
        
 if(CONC)      c.(caase){i} = h5read(currentfile,['/Conc/0']); 
        mass.(caase){i} = mean(mean(mean(c.(caase){i}))); 
        end;
        time.(caase){i}= h5read(currentfile,['/time'])
        
    end
    itime=i;
xu.(caase)=   h5read(currentfile,['/grid/xu']); 
xc.(caase)=   h5read(currentfile,['/grid/xc']);   
yv.(caase)=   h5read(currentfile,['/grid/yv']); 
yc.(caase)=   h5read(currentfile,['/grid/yc']); 
zw.(caase)=   h5read(currentfile,['/grid/zw']); 
zc.(caase)=   h5read(currentfile,['/grid/zc']); 



Lx=xu.(caase)(end);
Ly=yv.(caase)(end);
Lz=zw.(caase)(end);

timescale=pi*pi/Pe*((1/Lx).^2+(1/Ly).^2+(1/Lz).^2);

if(1== 1)
for l=1:length(time.(caase))
        
c_ana=[];
for i=1:length(xc.(caase)-1)
    for j=1:length(yc.(caase)-1)
        for k=1:length(zc.(caase)-1)
c_ana(i,j,k)=(1+cos(pi/Lx*xc.(caase)(i)).*cos(pi/Ly*yc.(caase)(j)).*cos(pi/Lz*zc.(caase)(k)).*exp(-timescale*time.(caase){l}))./2;
        end
    end
end
errorl2.(caase)(l)= sqrt( mean(mean(mean(((c.(caase){l}-c_ana).^2)  ))))

end

elseif(1==0 )
c_ana=[];
for i=1:length(xc.(caase)-1)
    for j=1:length(yc.(caase)-1)
        for k=1:length(zc.(caase)-1)
c_ana(i,j,k)=1+(sin(pi/Lx*xc.(caase)(i)).*sin(pi/Ly*yc.(caase)(j)).*sin(pi/Lz*zc.(caase)(k)).*exp(-timescale*time.(caase){itime}));
        end
    end
end
errorl2.(caase)= sqrt( mean(mean(mean(((c.(caase){itime}-c_ana).^2)  ))))


elseif(0)
    
   for l=1:length(time.(caase))
       
timescale=pi*pi/Pe*((2*2/Lx).^2);          
c_ana=[];
itime=l;

for i=1:length(xc.(caase)-1)
    for j=1:length(yc.(caase)-1)
        for k=1:length(zc.(caase)-1)
c_ana(i,j,k)=(sin(2*2*pi/Lx*(xc.(caase)(i)-0.2*time.(caase){itime})).*exp(-timescale*time.(caase){itime}));
        end
    end
end   

errorl2.(caase)(l)= sqrt( mean(mean(mean(((c.(caase){itime}-c_ana).^2)  ))))

   end
    
    
end

end



lineStyles = linspecer(6);

%%
timeindex=4

lineStyles = linspecer(6);
hfig=figure;
width = 2.5;     % Width in inches
height = 2.;     % Height in inches


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

tmp=c.(caase){timeindex}(:,25,25) ; 


plot( xc.(caase)  , tmp(:) ,  '-k','LineWidth', 1, 'MarkerSize', 3,'Color',lineStyles(1,:))
%hold on
%plot(  xc.(caase) , c_ana(:,1,1), '-k','LineWidth', 1, 'MarkerSize', 3,'Color',lineStyles(2,:) )



%%
x=xc.(caase)
y=yc.(caase)
z=zc.(caase)

yelo=yv.(caase)

figure
plot(y(:), pd(20,3:end-1, 20) )
hold on
plot(y(:), conc(20,3:end-1, 20), 'xr' )





figure
plot(yelo, Vunp(20,3:end-2, 20), 'ok' )
hold on
plot(yelo, Vp(20,3:end-2 ,20), 'xr' )



%%
x=xc.(caase)
y=yc.(caase)
z=zc.(caase)

yelo=yv.(caase)

figure
plot(y(:), pd(20,3:end-1, 20) )
hold on
plot(y(:), conc(20,3:end-1, 20), 'xr' )





figure
plot(yelo, Vunp(20,3:end-2, 20), 'ok' )
hold on
plot(yelo, Vp(20,3:end-2 ,20), 'xr' )

%%

lineStyles = linspecer(6);
hfig=figure;
width = 2.5;     % Width in inches
height = 2.;     % Height in inches


set(0,'defaulttextinterpreter', 'Latex')
fig = gcf;    
pos = get(gcf, 'Position');64
set(gcf, 'Position', [pos(1)*1.1 pos(2)*1.1 width*100, height*100]); %<- Set size
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/1.8;
bottom = (papersize(2)- height)/1.8;
myfiguresize = [left, bottom, width*1.1, height*1.1];
set(gcf,'PaperPosition', myfiguresize);

figure;
plot([time.(caase){:}], errorl2.(caase), 'x'  )


%%

timeindex=5
time.heat_exp{timeindex}

for i=1:length(xc.(caase)-1)
    for j=1:length(yc.(caase)-1)
        for k=1:length(zc.(caase)-1)
c_ana(i,j,k)=(sin(2*2*pi/Lx*(xc.(caase)(i)-0.2*time.(caase){timeindex})).*exp(-timescale*time.(caase){timeindex}));
        end
    end
end   

lineStyles = linspecer(6);
hfig=figure;
width = 2.5;     % Width in inches
height = 2.;     % Height in inches


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

tmp=c.(caase){timeindex}(:,5,5) ; 


plot( xc.(caase)  , tmp(:) ,  '-k','LineWidth', 1, 'MarkerSize', 3,'Color',lineStyles(1,:))
hold on
plot(  xc.(caase) , c_ana(:,1,1), '-k','LineWidth', 1, 'MarkerSize', 3,'Color',lineStyles(2,:) )

xlabel('$x$','Fontsize',9, 'Interpreter', 'Latex');
ylabel('$c$','Fontsize',9, 'Interpreter', 'Latex');
set(gca,'Fontsize',8,'YGrid', 'on', 'XGrid', 'on','XLim', [0 1]);%, 'XLim', [0 0.4], 'YLim', [0 400]);%, 'YScale', 'Log', 'YLim', [1e-12 1e-5]);%, 'YTick', [1e-10 1e-9 1e-8 1e-7 1e-5 1e-4]);
h=legend( 'parties', 'solution');
set(h,'Fontsize',8,'Location','NorthEast', 'Box', 'Off', 'Interpreter', 'none' );


print('-depsc2','-loose','-painters', [picturepath 'diffusiontest/convecsin_1']);

%%

timeindex=6;

tmp=u.(caase){timeindex}(:,1,1) ; 
hfig=figure;

plot( xc.(caase)  , tmp(:) ,  '-k','LineWidth', 1, 'MarkerSize', 3,'Color',lineStyles(1,:))













%%

timeindex=5;
%%%
hfig=figure;
width = 2.5;     % Width in inches
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


[Xcontour,Ycontour] = meshgrid(xc.(caase) , yc.(caase));


contourf( Xcontour, Ycontour,c.(caase){timeindex}(:,:,50)' ,150, 'LineStyle', 'None')
%hold on
%quiver( Xq(:) , Zq(:),  vx(:)*fac, vz(:)*fac, 0 ,'Color', [0 0 0], 'AutoScale', 'Off', 'LineWidth', 0.5)
%set(gca,'Fontsize',7, 'XLim', [0 0.75], 'YLim', [ -0.2 0.2], 'YTick', [-0.2 -0.1 0 0.1 0.2], 'XTick', [0 0.25 0.5 0.75]);
%quiver(0.7, 0.18 , 4*fac ,0, 0 ,'Color',  [0. 0. 0.], 'AutoScale', 'Off', 'LineWidth', 1.2,'MaxHeadSize', 1)

%annotation(hfig,'textbox',...
 %   [0.78 0.84 0.0442176870748304 0.0374305555555556],...
  %  'String','4',...
   % 'FitBoxToText','On',...
    %'EdgeColor','none', 'Color',  [0 0 0]);

xlabel('x','Fontsize',8);
ylabel('z','Fontsize',8, 'Interpreter','None');
%caxis([0 1.])
%num=30;
%tmp=hot(num)
%colormap(tmp(6:end,:) );

%colorbar('Fontsize', 8, 'YTick', [0 0.2 0.4 0.6 0.8 1])



%print(  [picturepath 'cosdiffusion/xx'],'-depsc','-painters', '-loose' )



%%

timeindex=6
lineStyles = linspecer(6);
hfig=figure;
width = 2.5;     % Width in inches
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

tmp=c.(caase){timeindex}(1,:,1) ; 

plot( yc.(caase)  , tmp(:) ,  '-k','LineWidth', 1, 'MarkerSize', 3,'Color',lineStyles(1,:))
hold on
plot(  yc.(caase) , c_ana(1,:,1), '-k','LineWidth', 1, 'MarkerSize', 3,'Color',lineStyles(2,:) )

xlabel('$x$','Fontsize',9, 'Interpreter', 'Latex');
ylabel('$y$','Fontsize',9, 'Interpreter', 'Latex');
set(gca,'Fontsize',8,'YGrid', 'on', 'XGrid', 'on');%, 'XLim', [0 0.4], 'YLim', [0 400]);%, 'YScale', 'Log', 'YLim', [1e-12 1e-5]);%, 'YTick', [1e-10 1e-9 1e-8 1e-7 1e-5 1e-4]);
h=legend( 'parties', 'solution');
set(h,'Fontsize',8,'Location','NorthEast', 'Box', 'Off', 'Interpreter', 'none' );
%matlabfrag( [picturepath 'shear/c1_poiseuille_uprof'], 'epspad', [0 ,10,0,0]);
%print('-depsc2','-loose','-painters', [picturepath 'shear/c1_poiseuille_uprof_c']);
%print(  [picturepath 'cosdiffusion/xx'],'-depsc','-painters', '-loose' )


