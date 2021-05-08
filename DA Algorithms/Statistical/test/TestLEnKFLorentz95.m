clear;
close all;

cd ../dev/

tic

% solution
K = 500;
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
J_H = eye(n);

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
y_sets = 2;
y = zeros(n, K-1, y_sets);

parfor i = 2:K+1
    for j=1:y_sets
        y(:, i-1, j) = H(i, x(:,i), R);
    end
end

H = @(t, x) x;
selector = cell(1,40);

parfor i1=1:n
    selector{i1} = simpleSelector(i1, n, 4);
end
sel = @(i) selector{i};

min = 41;
max = 50;

[ ARMSE_b, RMSE_b ] = averageRootMeanSquareError(x_b, x);
ARMSE_EnKF = zeros(y_sets, max);
ARMSE_EnKFLoc = zeros(y_sets, max);

ens = min:max;
for i=ens
    x_0b_en = ensembleInit(x_0b, P_0b, i);
    parfor j=1:y_sets
        y_cur = y(:,:,j);
        [ x_a_EnKF, x_a_EnKF_en] = da_seq_StochasticEnsembleKalmanFilter(x_0b_en, y_cur, @(t,x) M(t,x,F,T), H, R);
        [ x_a_EnKFLoc, x_a_EnKFLoc_en] = da_seq_StochasticEnsembleKalmanFilter(x_0b_en, y_cur, @(t,x) M(t,x,F,T), H, R, 1:K+1, 'Localization', sel);
        
        [ ARMSE_EnKF(j, i), RMSE_EnKF] = averageRootMeanSquareError(x_a_EnKF, x);
        [ ARMSE_EnKFLoc(j, i), RMSE_EnKFLoc] = averageRootMeanSquareError(x_a_EnKFLoc, x);
        j
    end
    i
    save Lorenz95Test
end

ARMSE_EnKF = mean(ARMSE_EnKF);
ARMSE_EnKFLoc = mean(ARMSE_EnKFLoc);

figure
plot(0, ARMSE_b, 'k*', ens, ARMSE_EnKF(ens), 'c*--', ens, ARMSE_EnKFLoc(ens), 'c*-');
legend('background ARMSE', 'EnKF ARMSE', 'EnKF ARMSE Localization');

time = toc

cd ../test
save Lorenz95Test