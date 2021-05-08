% Data Generation
clear
close all

burn_in = 5000;
K = burn_in + 10000;
T = 0.01;

radius = 1;
samples_circle = 100;
tetha_dot = 2*pi/(samples_circle*T); % After 100 T makes a complete circle

n = 2;

% Model
M = @(i, x) circle(i, x, tetha_dot, T);
J_M = @(i, x) finiteDifferenceAdjoint(i, x, M);

% Observer
J_H = @(t, x) eye(n);
H = @(t, x) eye(n);

x_0 = [ 10; 0 ];
x = zeros(n, K);
x(:, 1) = x_0;

for i=2:K
    x(:, i) = M(i, x(:, i-1));
end


R = 4*eye(n)
y = zeros(n, K-1);
for i=1:K-1
    y(:, i) = mvnrnd(x(:, i), R);
end

%
% Forget Everything
%
y_noise = zeros(size(y));
t = (1:K-1)';
for h=1:n        
    y_fit = fit(t, y(h, :)', 'smoothingspline', 'SmoothingParam', 0.1);
	y_noise(h, :) = y(h, :) - y_fit(t)';
end    

R = cov(y_noise');

uncertaintnyP = 2;
x_0b = mvnrnd(x_0, uncertaintnyP*R)';

H = @(t, x) x;

min = 2;
max = 25;

uncertaintnyP = 3;
x_0_en = cell(1,max);
parfor h=min:max
    x_0_en{h} = ensembleInit(x_0b, uncertaintnyP*R, h);
end

uncertaintnyP = 4;
P0b = uncertaintnyP*R;

clear min max
save TestData