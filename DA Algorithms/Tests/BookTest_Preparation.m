% Book Test:
% this test is like the one proposed throughtout the book.

% Preparation Script.

clear;
close all;

% solution
burn_in = 5*10^3;
K = burn_in + 10^4;

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
clear x_prev

t = 0:T:((K-1)*T);
t_inst = 2:K;

M = @(t, x) discreteLorenz95(t, x, F, T);
H = @(t, x, R) x + mvnrnd(zeros(length(x),1), R)';

%%% RUN the model, to get the true result (given the discretization).
r = 1;
P_0b = r*eye(n);
R = r*eye(n);    % Assume this small covariance matrix.

x = zeros(n, K);
x(:,1) = x0;
for i=2:K
    x(:,i) = M(i, x(:,i-1));
end

%%% MEASUREMENTS from the true results (given R).
y = zeros(n, K-1);

for i = 2:K
    y(:, i-1) = H(i, x(:,i), R);
end

%%% FORGET EVERYTHING: find x0, P0b, R
y_noise = zeros(size(y));
for h=1:n
    chunks = 1:100:K-1;
    chunks = [chunks K-1];
    for hh=2:length(chunks)
        tfit = (chunks(hh-1):chunks(hh))';
        y_fit = fit(tfit, y(h, tfit)', 'smoothingspline', 'SmoothingParam', 0.1);
        y_noise(h, tfit) = y(h, tfit) - y_fit(tfit)';        
    end
    h
end

R = cov(y_noise');

uncertaintnyP = 2;
x_0b = H(t, x0, uncertaintnyP*R);

H = @(t, x) x;

min = 2;
max = 250;

uncertaintnyP = 4;
x_0_en = cell(1,max);
parfor h=min:max
    x_0_en{h} = ensembleInit(x_0b, uncertaintnyP*R, h);
end

clear min max

save BookTestData