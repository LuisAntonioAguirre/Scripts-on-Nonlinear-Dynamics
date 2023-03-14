% Rossler discretized model

% LAA 03/08/21

clear 
close all

a = 0.398;
b = 2; 
c = 4; 

% initial conditions
x(1) = 0.1;
y(1) = 0.1;
z(1) = 0.1;

% discretization time
h = 0.05;

for k=1:10000
    x(k+1) = x(k)-h*(y(k)+z(k));
    y(k+1) = (1+a*h)*y(k) + h*x(k);
    z(k+1) = z(k) + h*(b+z(k)*(x(k)-c));
end

plot3(x,y,z)
