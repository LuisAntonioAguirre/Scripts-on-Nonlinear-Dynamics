using Plots
plotlyjs()

#=
  function dvLi(x,u)

implements the vector field for the Li system. First published in

Li, D., A three-scroll chaotic attractor, Physics Letters A, 372(4):387-393, 2008.

Default parameter values are a=42; c=11/6; d=0.16; e=0.65; k=55; F=20;


x state vector
u control input. Here it is placed in the first equation. This can be changed manually
xd time derivative of x (vector field at x). Because the original Li system is autonomous
we take u=0.

Luis A Aguirre 15/8/18
=#

function dvLi(x,u)
    a=42; c=11/6; d=0.16; e=0.65; k=55; F=20;
    xd=zeros(3)
    # Differential equations
    xd[1]=a*(x[2]-x[1])+d*x[1]*x[3]+u
    xd[2]=k*x[1]+F*x[2]-x[1]*x[3]
    xd[3]=c*x[3]+x[1]*x[2]-e*x[1]^2

    return xd
end


#=
   function rk(x0,u,h,t)

implements 4th-order Runge Kutta numerical integration algorithm
x0 is the state vector (before calling the integration function - it is the initial condition of the current step)
u is the control input, considered constant throughout the integration step h
h integration step (constant)
t is the time instant just before calling the integration function
=#

function rkLi(x0,u,h,t)
    # 1st call
    xd=dvLi(x0,u)
    savex0=x0
    phi=xd
    x0=savex0+0.5*h*xd

    # 2nd call
    xd=dvLi(x0,u)
    phi=phi+2*xd
    x0=savex0+0.5*h*xd

    # 3rd call
    xd=dvLi(x0,u)
    phi=phi+2*xd
    x0=savex0+h*xd

    # 4th call
    xd=dvLi(x0,u)
    x=savex0+(phi+xd)*h/6

    return x
end


# initial time
t0=0;
# final time
tf=50;
# integration interval
h=0.001;
# time vector
t=t0:h:tf;

# initial conditions
x0=[0.1, 0.1, 0.1]
# state evolution will be placed in x
x=[x0 zeros(length(x0),length(t)-1)]
# this system has no inputs
u=zeros(1,length(t))
#u=0.0*sin(1.01*t)

for k=2:length(t)
    x[:,k]=rkLi(x[:,k-1],u[k],h,t[k])
end

plot(x[1,:],x[2,:],x[3,:])

#writecsv("Li_data.txt",x[:,5000:end]')
