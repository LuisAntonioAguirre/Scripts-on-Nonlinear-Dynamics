
#=

Equation numbers in comments are from the paper:
Luis Antonio Aguirre and Christophe Letellier
Observability of multivariate differential embeddings
J. Phys. A: Math. Gen. 38 (2005) 6311-6326
doi:10.1088/0305-4470/38/28/004

see also
Relation between observability and differential embeddings for nonlinear dynamics
Christophe Letellier, Luis A. Aguirre, Jean Maquet
PHYSICAL REVIEW E 71, 066213 (2005)
doi: 10.1103/PhysRevE.71.066213

Luis A Aguirre 17/1/18
http://www.researcherid.com/rid/A-2737-2008

=#

function dvRossler(x,u)
    a=0.398; b=2; c=4 # spiral
#    a=0.52; # funnel
    xd=zeros(3)
    # Differential equations
    xd[1]=-x[2]-x[3]
    xd[2]=x[1]+a*x[2]+u
    xd[3]=b+x[3]*(x[1]-c)

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

function rkRossler(x0,u,h,t)
    # 1st call
    xd=dvRossler(x0,u)
    savex0=x0
    phi=xd
    x0=savex0+0.5*h*xd

    # 2nd call
    xd=dvRossler(x0,u)
    phi=phi+2*xd
    x0=savex0+0.5*h*xd

    # 3rd call
    xd=dvRossler(x0,u)
    phi=phi+2*xd
    x0=savex0+h*xd

    # 4th call
    xd=dvRossler(x0,u)
    x=savex0+(phi+xd)*h/6

    return x
end


# initial time
t0=0;
# final time
tf=500;
# integration interval
h=0.01;
# time vector
t=t0:h:tf;

# initial conditions
x0=randn(3,1)*0.1;
# state evolution will be placed in x
x=[x0 zeros(length(x0),length(t)-1)]
# this system has no inputs
u=zeros(1,length(t));
xd=zeros(3)

for k=2:length(t)
    # simulate the system
    x[:,k]=rkRossler(x[:,k-1],u[k],h,t[k]);
end
(linha,coluna)=size(x);


C1=[1 0 0];
C2=[0 1 0];
C3=[0 0 1];

# the following parameters must be the same as those in function dv_Rossler
a=0.389; b=2; c=4;


do1t=zeros(coluna,1)
do2t=zeros(coluna,1)
do3t=zeros(coluna,1)

# Evaluates \delta along a trajectory x(t): see Eq. (13)
# starting from N to avoid transients
N=2000;
	for i=N:coluna
        # jacobian matrix of Rossler system
        # Eq. (9) with n=0
		A=[0 -1 -1;1 a 0;x[3,i] 0 x[1,i]-c];

        # Definition using Lie Derivatives
        # Eq. (9) with n=1
        At2=[-1-x[3,i] -a c-x[1,i];a a^2-1 -1;b+2*x[3,i]*(x[1,i]-c) -x[3,i] (x[1,i]-c)^2-x[2,i]-2*x[3,i]];


		# observable x

        Qt1=[C1;C1*A;C1*At2]; # Eq. (8)
		Ot1=Qt1'*Qt1; # Argument in Eq. (13)
        (L1,v1)=eig(Ot1);
        vet1=sort(L1); # first eigenvalue is lambda_min and last will be lambda_max
		do1t[i]=abs(vet1[1])/abs(vet1[length(vet1)]); # Eq. (13)

		# observable y

        Qt2=[C2;C2*A;C2*At2]; # Eq. (8)
        Ot2=Qt2'*Qt2; # Argument in Eq. (13)
        (L2,v2)=eig(Ot2);
	    vet2=sort(L2); # first eigenvalue is lambda_min and last will be lambda_max
	    do2t[i]=abs(vet2[1])/abs(vet2[length(vet2)]); # Eq. (13)


		# observable z

        Qt3=[C3;C3*A;C3*At2]; # Eq. (8)
		Ot3=Qt3'*Qt3; # Argument in Eq. (13)
        (L3,v3)=eig(Ot3);
        vet3=sort(L3); # first eigenvalue is lambda_min and last will be lambda_max
		do3t[i]=abs(vet3[1])/abs(vet3[length(vet3)]); # Eq. (13)


	end; # i loop

    # average of doit along a trajectory (excluding possible initial
    # transients)

  	Do1t=mean(do1t[2000:coluna]) # Eq. (14)
  	Do2t=mean(do2t[2000:coluna]) # Eq. (14)
  	Do3t=mean(do3t[2000:coluna]) # Eq. (14)


    @printf "The observability coefficient for the variable x is: %s\r" round(Do1t,4)
    @printf "The observability coefficient for the variable y is: %s\r" round(Do2t,4)
    @printf "The observability coefficient for the variable z is: %s\r" round(Do3t,4)
