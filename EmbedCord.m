% Ilustration of choice of embedding dimension for the cord attractor

% LAA 8/02/21

clear; close all

%% 
% Simulates the cord system

t0=0; % initial time
tf=50; % final time
% integration interval
h=0.01;
t=t0:h:tf; % time vector


% initial conditions
x0=[0.1;0.1;0.1];

% initialization
x=[x0 zeros(length(x0),length(t)-1)];
% no inputs (no external force)
u=zeros(length(t),1);


for k=2:length(t)
    x(:,k)=rkCord(x(:,k-1),u(k),u(k),h,t(k));
end

figure(1)
plot3(x(1,1000:end),x(2,1000:end),x(3,1000:end),'k');
set(gca,'FontSize',18)
xlabel('x(t)')
ylabel('y(t)')
zlabel('z(t)')
axis([-6    6   -20    35  -25    35])
view([15.3 37.2])
grid


%%  
% Let us choose the time delay using mutual information
mutual(x(1,1000:end));
axis([-1 21 1 2.6])
delay=15;

%%
% Let us choose the dimension using the false nearest neighbors algorithm
Z=false_nearest(x(1,1000:end),1,5,delay);
false = Z(:,1:2);
figure(3)
plot(false(:,1),false(:,2),'ko-');
set(gca,'FontSize',18)
xlabel('dimension');
ylabel('FNN');

dim=3;

%%
% reconstruct the embedding space from a time series
y = phasespace(x(1,1000:end),dim,delay);

figure(4)
plot3(y(:,1),y(:,2),y(:,3),'k-','LineWidth',1);
set(gca,'FontSize',18)
xlabel('x(t)');
ylabel('x(t+\tau)');
zlabel('x(t+2\tau)');
axis([-5    10   -6  6   -5  6])
view([80.9 30.8])
grid


%%
%=========================================================
function xdot=dvCord(x,ux,uy,t)
%
% function xdot=dvCord(x,ux,uy,t)
% implements the Cord system
% x state vector
% the system has no inputs so ux=uy=0
% xd time derivative of x (vector field at x)

% The cord system

%
% First published in
%
% Aguirre, L.A., Letellier, C., Investigating observability properties from data in nonlinear dynamicsï¿œ, 
% Phisical Reivew E, 83:066209, 2011. DOI: 10.1103/PhysRevE.83.066209.
%
% Analysed in
%
% Letellier, C., Aguirre, L.A., Required criteria for recognizing new types of chaos: Application to the ï¿œcordï¿œ attractorï¿œ, 
% Phisical Reivew E, 85:036204, 2012. DOI: 10.1103/PhysRevE.85.036204.

a=0.258;
b=4.033;
F=8;
G=1;

xd(1)=-x(2)-x(3)-a*x(1)+a*F;
xd(2)=x(1)*x(2)-b*x(1)*x(3)-x(2)+G;
xd(3)=b*x(1)*x(2)+x(1)*x(3)-x(3);

xdot=xd';

% end dvCord
end


%=========================================================
function x=rkCord(x0,ux,uy,h,t)
% function x=rkCord(x0,ux,uy,h,t)
% 
% This function implements a numerical integration algorithm known as 4th-order Runge-Kutta  
% x0 is the state vector BEFORE calling the function (i.e. the initial condition at each integration step) 
% ux and uy (if different from zero) are external forces ("control actions") added to the first and second
% state equations (in dvExample.m), respectively. It is assumed that such
% control actions do NOT change during the integrarion period h. If
% ux=uy=0, the autonomous case is simulated.
% h integration interval.
% t the time BEFORE calling the function.
% The vector field (the equations) are in dvCord.m

% LAA 15/8/18

% 1st evaluation
xd=dvCord(x0,ux,uy,t);
savex0=x0;
phi=xd;
x0=savex0+0.5*h*xd;

% 2nd evaluation
xd=dvCord(x0,ux,uy,t+0.5*h);
phi=phi+2*xd;
x0=savex0+0.5*h*xd;

% 3rd evaluation
xd=dvCord(x0,ux,uy,t+0.5*h);
phi=phi+2*xd;
x0=savex0+h*xd;

% 4th evaluation
xd=dvCord(x0,ux,uy,t+h);
x=savex0+(phi+xd)*h/6;

% end rkCord
end
