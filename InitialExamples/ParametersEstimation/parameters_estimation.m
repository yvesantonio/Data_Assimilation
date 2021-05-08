clear;
close all;

%%% 
% 0. Build a test Model
%%%

u     = @(x) sin(pi*x);
du    = @(x) pi*cos(pi*x);
ddu   = @(x) -pi^2*sin(pi*x);

b_real  =   12;
c_real  =   0.4;

% IO
% -b*u'' + c*u' = f

%   b=12 , c=0.4
f = @(x) 12*pi^2*sin(pi*x) + 0.4*pi*cos(pi*x);

x   =	linspace(0, 1, 100);

% u_obs = u(t);
% Perturb u as if it was measured.
u_obs   =   u(x) + random('norm',0,0.01,1,length(x));

% u vs observations of u
figure;
plot(x,u(x), x,u_obs);

%%%
% 1.Solve the problem with a numerical solver
%%%
%%%
% 2. Solve for p
%%%
%%%
% 3. Try to solve the optimization problem with the steepest descent method
%%%

%
% 0. Initialize
% 1. In J.m
% 2. In DJ.m
% 3. use descentmethod.m
%

%%%
% Initialize
%%%
save('tmp.mat','u_obs','f');

b0  = 5;
c0  = 1;

f_J = @(m) J(m, x);

J([12,0.4], x)
DJ([12,0.4])

J([10,2], x)
DJ([10,2])

[mk, Jk, counter, error, mks, Jks] = descentmethod(f_J,[b0;c0],1e-5,100);

figure
plot(mks(:,1),mks(:,2), mk(1),mk(2),'*')

J(mk, x)
DJ(mk)