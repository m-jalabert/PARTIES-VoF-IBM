clear;
clc,

% Geometry parameters
h = 1;

% Velocity of plate
u_0 = 8; % define this better

y = (0:0.01:1)'; % column vector
velocity=zeros(101,1); % column vector

for i = 1:101
    velocity(i) = u_0 * (y(i)/h);
end

% File handling
fileID = fopen('analytical_velo_CF.dat','w');
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
legend('U(y)')
title('Analytical Solution')
