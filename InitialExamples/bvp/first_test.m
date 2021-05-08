clear;
close all;

A = [   0   1;
        0   0.4/12; ];
    
B = [   0;  -1/12];

f = @(t) 12*pi^2*sin(pi*t) + 0.4*pi*cos(pi*t);


%%% Initial conditions ????
%%% ?? The book uses boundary conditions on u, not initial conditions.

odefun  = @(t, x) A*x + B*f(t);
res     = @(in, fin) [in(1) fin(1)];
solinit = bvpinit(linspace(0, 1, 10),[0 0]);

sol  = bvp4c(odefun, res, solinit);
