using Plots
plotlyjs()

#=
  function dvKuramoto(x,n,w,rho,Adj,t)

  implements the equations for a network of Kuramoto oscillators
  x is the state of the full network of size n
  w is the (column) vector of parameters (natural frequency of each oscillator)
  the number of nodes in the network is length(w)=length(x)=n
  rho is the coupling strength which is the same for every existing
  coupling
  Adj is the adjacency matrix. It indicates the connections among the
  oscillators. Adj(i,j)=1 indicates that there is a directed link from node
  j to node i. This is a square matrix of dimension length(x)=n.
  t is the time vector
  xd is the time derivative of x (vector field at x)

  LAA 8/1/19
=#

function dvKuramoto(x,n,w,rho,Adj,t)

    xd=w + rho*sum(sin(ones(n,1)*x' -x*ones(1,n)).*Adj, 2);

    return xd
end


#=
    function x=rkKuramoto(x0,n,w,rho,Adj,h,t)

    implements numerical integration using 4th-order Runge Kuta
    x0 state vector before call the function (initial conditions)
    h integration interval
    t time instant just before calling the function
    the vector field is in dvKuramoto.m

=#

function rkKuramoto(x0,n,w,rho,Adj,h,t)
    # 1st call
    xd=dvKuramoto(x0,n,w,rho,Adj,t)
    savex0=x0
    phi=xd
    x0=savex0+0.5*h*xd

    # 2nd call
    xd=dvKuramoto(x0,n,w,rho,Adj,t+0.5*h)
    phi=phi+2*xd
    x0=savex0+0.5*h*xd

    # 3rd call
    xd=dvKuramoto(x0,n,w,rho,Adj,t+0.5*h)
    phi=phi+2*xd
    x0=savex0+h*xd

    # 4th call
    xd=dvKuramoto(x0,n,w,rho,Adj,t+h)
    x=savex0+(phi+xd)*h/6

    return x
end


# initial time
t0=0;
# final time
tf=200;
# integration interval
h=0.1;
# time vector
t=t0:h:tf;


# Graph of the network
n = 10;                                                      # number of nodes
Adj = diagm(ones(Int8,n-1),-1) + diagm(ones(Int8,n-1),+1);  # adjacency matrix
rho = 0.15;                                                  # coupling

# initial conditions
x0=rand(n)
# natural frequencies of the oscillators taken randomly around 1 rad/s
w=randn(n)*0.01+ones(n)
# state evolution will be placed in x
x=[x0 zeros(n,length(t)-1)]
xd=zeros(n)

for k=2:length(t)
    x[:,k]=rkKuramoto(x[:,k-1],n,w,rho,Adj,h,t[k])
end

# vector of absolute phase differences between oscillators o1 and o2
o1=1;
o2=n;
errp=abs(x[o1,:]-x[o2,:])

plot(t,x[1,:],yaxis="x_1, x_n", xaxis="time")
plot!(t,x[n,:])

plot(t,errp,yaxis="|x_1-x_n|",xaxis="time")
