function [dJ] = DJ(m)

% J and gradient J
% J     = 1/2*integral[0,1]((u-u_obs)^2*dx)
% dJ    = integral[0,1]((u-u_obs)* u^ * dx)
% u^    = (u~ - u)/alpha;       alpha -> 0
% AFTER COMPUTATIONS:
% dJ    = [ integral[0,1](p*ddu*dx)
%           integral[0,1](-p*du*dx) ];

% Adjoint problem
% IO
% -b*p'' - c*p' = (u-u_obs)


%%%
% 2. Solve for p
%%%

load tmp

sol  = solve_direct([m(1),-m(2)], err, x_ode_1);

x_p = sol.x;
p   = sol.y(1,:);

%%%
% figure
% plot(t_p, p);
%%%

if(length(du) == length(p))
    f_b =   @(x) interp1(x_ode_1, p.*ddu, x);
    f_c =   @(x) interp1(x_ode_1, p.*du, x);
else
    warning('solver changed grid, again')
    p   =   @(x) interp1(sol.x, p, x);
    du  =   @(x) interp1(x_ode_1, du, x);
    ddu =   @(x) interp1(x_ode_1, ddu, x);
    
    f_b =   @(x) p(x).*ddu(x);
    f_c =   @(x) p(x).*du(x);
end
dJb = integral(f_b, 0, 1);
dJc = -integral(f_c, 0, 1);

dJ = [ dJb;
       dJc; ];