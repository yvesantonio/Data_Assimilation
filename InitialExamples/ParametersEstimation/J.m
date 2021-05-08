function J = J(m, x_ode)
% J     = 1/2*integral[0,1]((u-u_obs)^2*dx)
load tmp;

[sol] = solve_direct(m, f, x_ode);

u       =   sol.y(1,:);
du      =   sol.y(2,:);
ddu     =   sol.yp(2,:);

%%%
%   Is it right to use interp1?
%%%
x_ode_1 = x_ode;
if(length(u) == length(u_obs))
    err     =	@(x) interp1(x_ode, (u-u_obs), x);
else
    warning('solver changed grid')
    u_obs_1 =   @(x) interp1(x_ode, u_obs, x);
    x_ode_1 =	sol.x; 
    u       =	@(x) interp1(x_ode_1, u, x);
    err     =   @(x) u(x) - u_obs_1(x);
end
integrand = @(x) err(x).^2;
J   =	1/2*integral(integrand,0,1);

save('tmp.mat', 'f','u_obs', 'u', 'du', 'ddu', 'x_ode_1', 'err');