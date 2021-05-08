% Lorenz 95 Equations Script
clear; close all;

% solution
m = 500;
n = 40;
T = 0.05;

F = 8;
%%%
%   I don't want to struggle to get some initial conditions...
%%%
x0 = zeros(n,1);
%%%
%   I perturb the first term...
%%%
x0(1) = F/10;
[t, x_prev] = ode45(@(t,x) Lorenz95(t,x,F), [0,F*10], x0);

%%%
%   Here we are.
%%%
x0 = x_prev(end,:)';

t = 0:T:(m*T);
t_inst = 0:m;

[t, x] = ode45(@(t,x) Lorenz95(t,x,F), t, x0);
x = x';
x = [   x(end,:);
        x;];
    
[T_grid, X_grid] = meshgrid(t_inst,0:n);

figure
surf(T_grid,X_grid, x, 'edgecolor','none')
c = gray(256);
c = flipud(c);
colormap(c)
shading interp
colorbar
view(0,90)

%%%
%   Test the discrete model.
%%%

x_discrete = zeros(40, m+1);
x_discrete(:,1) = x0;
for i=2:m+1
    x_discrete(:,i) = discreteLorenz95(i, x_discrete(:,i-1), F, T);
end

x_discrete = [  x_discrete(end,:);
                x_discrete;];

figure
surf(T_grid,X_grid, x_discrete, 'edgecolor','none')
c = gray(256);
c = flipud(c);
colormap(c)
shading interp
colorbar
view(0,90)