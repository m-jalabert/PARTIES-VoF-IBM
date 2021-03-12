clear;
clc,

% Geometry parameters
h = 1;

% kinematic viscosity
nue = 1; 

% Pressure gradient
delta_p = -12.037; % define this better

y = linspace(0,1,100)'; % column vector
velocity=zeros(100,1); % column vector

% Calculate the velocity
for i = 1:100
    velocity(i) = -( (delta_p)/(2*nue) ) * (h - y(i)) * y(i);
end

% File handling
fileID = fopen('analytical_velo_PF.dat','w');
fprintf(fileID,'%10s\n','U(y)');
fprintf(fileID,'%10.8f\n',velocity);
fclose(fileID);

% plotting
y_axis = linspace(1,100,100);
plot(velocity,y_axis) 
grid on
grid minor
legend('U(y)')
title('Analytical Solution')
