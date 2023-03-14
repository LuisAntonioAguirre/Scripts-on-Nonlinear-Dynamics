% investigates the effect of a small parameter that
% breaks the symmetry of the normal equation for the
% pitchfork bifurcation. This should help you see that if symmetry is
% broken, instead of a pitchfork bifurcation we get: 1) a stable fixed 
% point that does not bifurcate plus the creation of a stable and an
% unstable fixed points via saddle-node bifurcation.

% LAA 10/10/2022

% this is the symmetry-breaking parameter
delta = 0.00001;

% bifurcation parameter
mu = -0.01:0.0001:0.05;

% location of fixed points
LF = zeros(3,length(mu));

for k = 1:length(mu)
    LF(:,k) = sort(real(roots([-1 0 mu(k) delta])));
end

figure(1)
plot(mu,LF,'k')
set(gca,'FontSize',16)
xlabel('\mu')
ylabel('x')

%% stability

clear mu
fp = zeros(3,1);

% value of mu of interest
mu = 0.01;
% fixed point at that value of mu
fp = roots([-1 0 mu delta])
% chose the fixed point whose stability is to be investigated
x = fp(2);

dFdx = -3*x^2 + mu;

if dFdx < 0
    disp('fixed point is stable');
else 
    disp('fixed point is unstable');
end



