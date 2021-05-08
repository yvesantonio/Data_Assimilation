function [ sol ] = solve_direct( m, f, x_ode)
% SS
% u''   =  c/b*u' - 1/b*f;
% u'    =  u';

A = [   0   1;
        0   m(2)/m(1); ];
    
B = [   0;  -1/m(1)];

odefun  = @(x, u) A*u + B*f(x);
res     = @(in, fin) [in(1); fin(1)];
solinit = bvpinit(x_ode,[0 0]);

sol  = bvp4c(odefun, res, solinit);

end