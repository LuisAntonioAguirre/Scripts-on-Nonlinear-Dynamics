using Plots
plotlyjs()

#=
  function dvCord(x,u)

implements the vector field for Cord system. First published in

Aguirre, L.A., Letellier, C., Investigating observability properties from data in nonlinear dynamics,
Phisical Reivew E, 83:066209, 2011. DOI: 10.1103/PhysRevE.83.066209.

Analysed in

Letellier, C., Aguirre, L.A., Required criteria for recognizing new types of chaos: Application to the cord attractor,
Phisical Reivew E, 85:036204, 2012. DOI: 10.1103/PhysRevE.85.036204.

Default parameter values are  a=0.258; b=4.033; F=8; G=1;
x state vector
u control input. Here it is placed in the first equation. This can be changed manually
xd time derivative of x (vector field at x)

Luis A Aguirre 11/12/17
http://www.researcherid.com/rid/A-2737-2008
=#

function dvCord(x,u)
    a=0.258; b=4.033; F=8; G=1;
    xd=zeros(3)
    # Differential equations
    xd[1]=-x[2]-x[3]-a*x[1]+a*F+u
    xd[2]=x[1]*x[2]-b*x[1]*x[3]-x[2]+G
    xd[3]=b*x[1]*x[2]+x[1]*x[3]-x[3]

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

function rkCord(x0,u,h,t)
    # 1st call
    xd=dvCord(x0,u)
    savex0=x0
    phi=xd
    x0=savex0+0.5*h*xd

    # 2nd call
    xd=dvCord(x0,u)
    phi=phi+2*xd
    x0=savex0+0.5*h*xd

    # 3rd call
    xd=dvCord(x0,u)
    phi=phi+2*xd
    x0=savex0+h*xd

    # 4th call
    xd=dvCord(x0,u)
    x=savex0+(phi+xd)*h/6

    return x
end


# initial time
t0=0;
# final time
tf=200;
# integration interval
h=0.01;
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
    x[:,k]=rkCord(x[:,k-1],u[k],h,t[k])
end

# plot(x[1,:],x[2,:],x[3,:])

# writecsv("cord_data.txt",x[:,2000:end]')
