% Lorenz discretized model

% LAA 03/08/21

clear 
close all

sigma = 10;
r = 28; 
b = 8/3; 

% initial conditions
x(1) = 9;
y(1) = -0.4;
z(1) = 36.5;

% discretization time
h = 0.01;

for k=1:5000
    x(k+1) = (1-sigma*h)*x(k)+h*sigma*y(k);
    y(k+1) = (1-h)*y(k) + h*x(k)*(r-z(k));
    z(k+1) = (1-b*h)*z(k) + h*x(k)*y(k);
end

plot3(x,y,z)