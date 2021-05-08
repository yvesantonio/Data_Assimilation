clear;
close all;

u     = @(x) sin(pi*x);
du    = @(x) pi*cos(pi*x);
ddu   = @(x) -pi^2*sin(pi*x);

b_real  =   12;
c_real  =   0.4;

% IO
% -b*u'' + c*u' = f

%   b=12 , c=0.4
f = @(x) 12*pi^2*sin(pi*x) + 0.4*pi*cos(pi*x);

x   = 0:0.01:1;
sol = solve_direct([12, 0.4], f, x);
x_ode       =   sol.x;
u_calc      =   sol.y(1,:);
du_calc     =   sol.y(2,:);
ddu_calc    =   sol.yp(2,:);

figure;
plot(x_ode, u(x_ode), x_ode, u_calc,'*');
figure
plot(x_ode, du(x_ode), x_ode, du_calc,'*');
figure
plot(x_ode, ddu(x_ode), x_ode, ddu_calc,'*');