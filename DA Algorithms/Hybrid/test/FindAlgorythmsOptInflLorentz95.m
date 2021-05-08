clear;
close all;

cd ../dev

tic
%% parameters

% solution
K = 400;
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

t = 0:T:(K*T);
t_inst = 0:K;

M = @(t, x, F, T) discreteLorenz95(t, x, F, T);
H = @(t, x, R) x + mvnrnd(zeros(length(x),1), R)';
%%% RUN the model, to get the true result (given the discretization).
%   Use the T = 0.05;

r = 4;
P_0b = r*eye(n);
R = r*eye(n);    % Assume this small covariance matrix.
Q = zeros(n);

x_0b = H(0, x0, R);
x_b = x_0b.*eye(n)*ones(n, K+1);
x = zeros(n, K+1);
x(:,1) = x0;

for i = 2:K
    x(:,i) = M(i, x(:,i-1), F, T);
    x_b(:,i) =	M(i, x_b(:,i-1), F, T);
end

%%% Test 1: assimilation every step.
y_sets = 5;
y = zeros(n, K-1, y_sets);

parfor i = 2:K+1
    for j=1:y_sets
        y(:, i-1, j) = H(i, x(:,i), R);
    end
end

% infl = 1:0.1:1.5;
% infl = 1:0.05:1.3;
% infl = 1:0.025:1.1;
infl = 1.015:0.01:1.045;


l_min = 27;
l_max = 30;

max = length(infl);
ARMSE_MLETKFInfl = zeros(l_max, y_sets, max);

for l=l_min:l_max
    x_0b_en = ensembleInit(x_0b, P_0b, l);
    for i=1:max
        parfor j=1:y_sets
            [ x_a_MLETKFInfl, x_a_MLETKFInfl_en] = da_hyb_MaximumLikelihoodEnsembleTransformKalmanFilter(x_0b_en, y(:,:,j), @(t,x) M(t,x,F,T), @(t,x) H(t, x, 0*eye(n)), R, 1:K+1, 'MaxIter', 25, 'Inflation', infl(i));
            [ ARMSE_MLETKFInfl(l, j, i), RMSE] = averageRootMeanSquareError(x_a_MLETKFInfl, x);
            RMSE_MLETKFInfl(l, j, i) = RMSE(end);
        end     
    end
    l
    save OptInfl95
end

ARMSE_MLETKFInfl_mean = squeeze(mean(ARMSE_MLETKFInfl, 2))';

[ens_grid, infl_grid] = meshgrid(l_min:l_max,infl);

figure
surf(ens_grid, infl_grid, real(ARMSE_MLETKFInfl_mean(:, l_min:l_max)))

xlabel('ensemble')
ylabel('infl')
zlabel('ARMSE')
%%%%
% BEST: ~1.025 Inflation (not so sure, run with an higher y_sets.
% Ensemble dimension: the higher the better.
%%%%


cd ../test
save OptInfl95