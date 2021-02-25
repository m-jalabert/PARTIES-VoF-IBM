clear;
clc,

% Geometry parameters
h = 1;

% Velocity of plate
u_0 = 8; % define this better

y = linspace(0,1,100)'; % column vector
velocity=zeros(100,1); % column vector

for i = 1:100
    velocity(i) = u_0 * (y(i)/h);
end

% File handling
fileID = fopen('analytical_velo.dat','w');
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
