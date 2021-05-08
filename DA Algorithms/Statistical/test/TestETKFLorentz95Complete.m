clear;
close all;

cd ../dev/

% solution
burn_in = 5*10^3;
K = burn_in + 10^5;

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

t = 0:T:((K-1)*T);
t_inst = 1:K;

M = @(t, x) discreteLorenz95(t, x, F, T);
H = @(t, x, R) x + mvnrnd(zeros(length(x),1), R)';

%%% RUN the model, to get the true result (given the discretization).
r = 1;
P_0b = r*eye(n);
R = r*eye(n);    % Assume this small covariance matrix.

x_0b = H(t,x0,P_0b);
x = zeros(n, K);
x(:,1) = x0;
for i=2:K
    x(:,i) = M(i, x(:,i-1));
end

%%% Test 1: assimilation every step.
y = zeros(n, K);

for i = 2:K
    y(:, i-1) = H(i, x(:,i), R);
end

min = 5;
max = 8;
ens = min:max;

H = @(t, x) x;
selector = cell(1,40);

parfor i1=1:n
    selector{i1} = roundSelector(i1, n, 4);
end
sel = @(i) selector{i};

ARMSE_ETKF = zeros(1, max);
ARMSE_ETKFInfl = zeros(1, max);
ARMSE_ETKFLoc = zeros(1, max);
ARMSE_ETKFInflLoc = zeros(1, max);
RMSE_ETKF = zeros(max, K - burn_in);
RMSE_ETKFInfl = zeros(max, K - burn_in);
RMSE_ETKFLoc = zeros(max, K - burn_in);
RMSE_ETKFInflLoc = zeros(max, K - burn_in);

x_0_en = cell(1,max);
parfor h=ens
    x_0_en{h} = ensembleInit(x_0b, P_0b, h);
end

'Initialized...'
tic
parfor i=ens
    [ x_a_ETKF, x_a_ETKF_en] = da_seq_EnsembleTransformKalmanFilter(x_0_en{i}, y, M, H, R);
    [ ARMSE_ETKF(i), RMSE_ETKF(i, :) ] = averageRootMeanSquareError(x_a_ETKF(:, burn_in+1:K), x(:, burn_in+1:K));
    
    [ x_a_ETKFInfl, x_a_ETKFInfl_en] = da_seq_EnsembleTransformKalmanFilter(x_0_en{i}, y, M, H, R, t_inst, 'Inflation', 1.025);
    [ ARMSE_ETKFInfl(i), RMSE_ETKFInfl(i, :)] = averageRootMeanSquareError(x_a_ETKFInfl(:, burn_in+1:K), x(:, burn_in+1:K));
    
    [ x_a_ETKFLoc, x_a_ETKFLoc_en] = da_seq_EnsembleTransformKalmanFilter(x_0_en{i}, y, M, H, R, t_inst, 'Localization', sel);
    [ ARMSE_ETKFLoc(i), RMSE_ETKFLoc(i, :)] = averageRootMeanSquareError(x_a_ETKFLoc(:, burn_in+1:K), x(:, burn_in+1:K));
    
    [ x_a_ETKFInflLoc, x_a_ETKFInflLoc_en] = da_seq_EnsembleTransformKalmanFilter(x_0_en{i}, y, M, H, R, t_inst, 'Inflation', 1.025, 'Localization', sel);
    [ ARMSE_ETKFInflLoc(i), RMSE_ETKFInflLoc(i,:)] = averageRootMeanSquareError(x_a_ETKFInflLoc(:, burn_in+1:K), x(:, burn_in+1:K));
end
time = toc

figure
plot(ens, ARMSE_ETKF(ens), 'g*--', ens, ARMSE_ETKFInfl(ens), 'g*-', ens, ARMSE_ETKFLoc(ens), 'c*-', ens, ARMSE_ETKFInflLoc(ens), 'c*--');
title('ARMSE');
legend('ETKF', 'ETKF Inflation', 'ETKF Localization', 'ETKF full');

cd ../test