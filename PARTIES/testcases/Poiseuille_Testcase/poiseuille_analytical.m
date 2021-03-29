clear;
clc,

% Geometry parameters
h = 1;

% kinematic viscosity
nue = 1; 

% Pressure gradient
delta_p = -12.037; % define this better

y = (0:0.01:1)'; % column vector
velocity=zeros(101,1); % column vector

% Calculate the velocity
for i = 1:101
    velocity(i) = -( (delta_p)/(2*nue) ) * (h - y(i)) * y(i);
end

% File handling
fileID = fopen('analytical_velo_PF.dat','w');
% fprintf(fileID,'%10s\n','U(y)');
% fprintf(fileID,'%10.8f\n',velocity);
fprintf(fileID,'%5s\t%10s\n','y','U(y)');
for i = 1:101
    fprintf(fileID,'%5.4f\t%10.8f\n',y(i),velocity(i));
end
fclose(fileID);

% plotting
y_axis = linspace(0,100,101);
plot(velocity,y_axis)
grid on
grid minor
legend('U(x)')
title('Analytical Solution')
