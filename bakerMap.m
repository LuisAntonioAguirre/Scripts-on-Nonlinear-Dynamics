% implements the baker's map

% LAA 27/09/18

close all
clear

% stretching and compression factor (default sc=2)
sc=2;

% number of iterates is N-1 (because the initial condition
% corresponds to N=1).
N=8;

% number of realizations. I you just want to iterate the map for one
% initial condition, choose Nr=1.
Nr=3000;

% Definition of matrices with the states for all iterations and realization.
x=zeros(N,Nr);
y=zeros(N,Nr);

% Nr initial values taken randomly from a uniform
% distribution within the unit square
x(1,:)=rand(1,Nr)/2;
y(1,:)=rand(1,Nr);


figure(1)
set(gca,'FontSize',18)
plot(x(1,:),y(1,:),'k.')
xlabel('x')
ylabel('y')
axis([0 1 0 1])


% k is the iterate
% i is de realization
for k=2:N
    for i=1:Nr
        if x(k-1,i)<1/2
            x(k,i)=sc*x(k-1,i);
            y(k,i)=y(k-1,i)/sc;
        else
            x(k,i)=sc*x(k-1,i)-1;
            y(k,i)=(y(k-1,i)+1)/sc;
        end
    end
    
    figure(k)
    set(gca,'FontSize',18)
    plot(x(k,:),y(k,:),'b.')
    xlabel('x')
    ylabel('y')
    axis([0 1 0 1])

end


%% selection of figures

 figure(k+1)
 subplot(221)
 set(gca,'FontSize',18)
 plot(x(1,:),y(1,:),'k.')
 title('(a)')
 xlabel('x')
 ylabel('y')
 axis('equal')
 axis([-0.05 1.05 -0.05 1.05])
 
 
 subplot(222)
 set(gca,'FontSize',18)
 plot(x(2,:),y(2,:),'k.')
 title('(b)')
 xlabel('x')
 ylabel('y')
 axis('equal')
 axis([-0.05 1.05 -0.05 1.05])

 subplot(223)
 set(gca,'FontSize',18)
 plot(x(4,:),y(4,:),'k.')
 title('(c)')
 xlabel('x')
 ylabel('y')
 axis('equal')
 axis([-0.05 1.05 -0.05 1.05])
            
 subplot(224)
 set(gca,'FontSize',18)
 plot(x(8,:),y(8,:),'k.')
 title('(d)')
 xlabel('x')
 ylabel('y')
 axis('equal')
 axis([-0.05 1.05 -0.05 1.05])
    
% print -dpng bakersMap.png