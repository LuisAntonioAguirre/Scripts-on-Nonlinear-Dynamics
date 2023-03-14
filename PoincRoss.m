% Poincare Section for Rossler system
% Luis Aguirre 05/10/17

clear
close all

% time definitions
t0=0;
tf=1500;
h=0.01;
t=t0:h:tf;

% if parameters should be changed, remember to change also in function rkRossler (see below)
a=0.398; b=2;c=4; % spiral
%a=0.52; % funnel type
e=randn(3,1)*0.1;
    
% Poincare section
% x at fixed point and dot(x)>0, this will be tested later
xp=(c-sqrt(c^2-4*a*b))/2;

% initial conditions
x0=[e(1); e(2); e(3)];

% initialization
x=[x0 zeros(length(x0),length(t)-1)];
% no inputs (no external force)
u=zeros(length(t),1);


% initial time do avoid transients
ti=(tf/2)/h;


% initialiation of vector with indices of values at Poincare section
kp=zeros(1000,1);
i=1;
for k=2:length(t)
    % simulates the system
    x(:,k)=rkRossler(x(:,k-1),u(k),h,t(k));       
end

k=ti;
i=1;
while k<length(t)
if x(1,k)>xp % if x passed fixed point coordinate 
        if -x(2,k)-x(3,k)>0 % and if dot(x)>0
            kp(i)=k;
            i=i+1;
            k=k+500;
        end
end  
k=k+1;
end





% without transients
X=x(:,ti:length(t));

figure(1)
plot3(X(1,:),X(2,:),X(3,:),'k');
hold on
plot3(x(1,kp(2:i-1)),x(2,kp(2:i-1)),x(3,kp(2:i-1)),'r.');
hold off
%plot(x(:,1),x(:,2));

%%

% interpolation

for j=2:i-1
    % inclination of interpolation line in x-y plane
    alphay=(x(2,kp(j))-x(2,kp(j)-1))/(x(1,kp(j))-x(1,kp(j)-1));
    % intercept of interpolation line in x-y plane
    betay=(x(2,kp(j)-1)*x(1,kp(j))-x(2,kp(j))*x(1,kp(j)-1))/(x(1,kp(j))-x(1,kp(j)-1));
    % line
    ypo(j)=betay+alphay*xp;
    
    % inclination of interpolation line in z-y plane
    alphaz=(x(3,kp(j))-x(3,kp(j)-1))/(x(2,kp(j))-x(2,kp(j)-1));
    % intercept of interpolation line in z-y plane
    betaz=(x(3,kp(j)-1)*x(2,kp(j))-x(3,kp(j))*x(2,kp(j)-1))/(x(2,kp(j))-x(2,kp(j)-1));
    % line
    zpo(j)=betaz+alphaz*ypo(j);
end

% original data and interpolated data
%[x(1,kp(1:i-1))' x(2,kp(1:i-1))' x(3,kp(1:i-1))' xp*ones(i-1,1) ypo' zpo']

figure(1)
set(gca,'FontSize',18)
plot3(X(1,:),X(2,:),X(3,:),'k');
xlabel('x');
ylabel('y');
zlabel('z');
grid
hold on
plot3(x(1,kp(2:i-1)),x(2,kp(2:i-1)),x(3,kp(2:i-1)),'ro',xp*ones(i-2,1),ypo(2:i-1)',zpo(2:i-1)','bv');
hold off

figure(2)
set(gca,'FontSize',18)
plot(X(1,:),X(2,:),'k');
xlabel('x(t)')
ylabel('y(t)')
hold on
plot(x(1,kp(2:i-1)),x(2,kp(2:i-1)),'ko',xp*ones(i-2,1),ypo(2:i-1)','bv',x(1,kp(2:i-1)-1),x(2,kp(2:i-1)-1),'ko');
hold off


%%

% Poincare section
figure(3)
set(gca,'FontSize',18)
plot(ypo(2:i-1),zpo(2:i-1),'k.')
xlabel('y_n') 
ylabel('z_n') 



% first-return map
figure(4)
set(gca,'FontSize',18)
plot(-ypo(2:i-2),-ypo(3:i-1),'k.',[0 12],[0 12],'k-')
xlabel('-y_n') 
ylabel('-y_{n+1}') 


%=======================================================
function xdot=dvRossler(x,u,t)
%
% function xd=dvRossler(x,u) 
% implements the vector field for Rossler's system
% x state vector
% u control input. Here it is placed in the second equation. This can be
% changed manually
% xd time derivative of x (vector field at x)

% Luis A Aguirre 3/6/16
% http://www.researcherid.com/rid/A-2737-2008

%a=0.398; b=2; c=4; %standard
a=0.2; b=0.2; c=5.7; %Pyragas paper

% Differential equations
xd(1)=-x(2)-x(3);
xd(2)=x(1)+a*x(2)+u;
xd(3)=b+x(3)*(x(1)-c);

xdot=xd';
end;



%=======================================================
function x=rkRossler(x0,u,h,t)
% function x=rk(x0,u,h,t)
% 
% implements 4th-order Runge Kutta numerical integration algorithm 
% x0 is the state vector (before calling the integration function - it is the initial condition of the current step)
% u is the control input, considered constant throughout the integration step h 
% h integration step (constant)
% t is the time instant just before calling the integration function 
% the vector field in in dvRossler.m

% Luis A Aguirre 3/6/16
% http://www.researcherid.com/rid/A-2737-2008

% 1st call
xd=dvRossler(x0,u,t);
savex0=x0;
phi=xd;
x0=savex0+0.5*h*xd;

% 2nd call
xd=dvRossler(x0,u,t+0.5*h);
phi=phi+2*xd;
x0=savex0+0.5*h*xd;

% 3rd call
xd=dvRossler(x0,u,t+0.5*h);
phi=phi+2*xd;
x0=savex0+h*xd;

% 4th call
xd=dvRossler(x0,u,t+h);
x=savex0+(phi+xd)*h/6;
end;








