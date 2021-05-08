clear;
close all;

x_b = 12;
x   = 10;
H   = 1;
B   = 5; %very high
R   = 0.1;
y   = H*x + random('normal', 0, R);

x_a_formula = x_b + B*H'*(H*B*H' + R)^-1*(y-H*x_b)

%   The book always suggest quasi-Newton methods.
%   Supply the gradient.
%   Supply max iterations
%   Supply toleranceù
solver_options = optimoptions(@fminunc, 'Algorithm', 'quasi-newton', 'SpecifyObjectiveGradient', true, 'MaxIterations', 100, 'OptimalityTolerance', 1e-6, 'Display', 'iter' );

cd ../dev
[ x_a, J_xa] = da_var_3DVar(x_b, y, B, R, H, @(f, x0) fminunc(f, x0, solver_options));


x_b = [ 12, 175 ]';
x   = [ 10, 142 ]';
H   = [ 1	1];
B   = [ 5	0;
        0   50; ];    
R   = 11;
y   = H*x + random('normal', 0, R);

x_a_formula = x_b + B*H'*(H*B*H' + R)^-1*(y - H*x_b)

[ x_a, J_xa] = da_var_3DVar(x_b, y, B, R, H, @(f, x0) fminunc(f, x0, solver_options))
