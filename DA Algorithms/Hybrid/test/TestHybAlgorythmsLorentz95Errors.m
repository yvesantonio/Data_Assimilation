clear;
close all;

cd ../dev

tic

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
x_dis = [   x(end,:);
            x;];

x_b_dis = [ x_b(end,:);
            x_b;];

[T_grid, X_grid] = meshgrid(t_inst,0:n);

figure
surf(T_grid,X_grid, x_dis, 'edgecolor','none')
legend('real')
c = gray(256);
c = flipud(c);
colormap(c)
shading interp
colorbar
view(0,90)

figure
surf(T_grid,X_grid, x_b_dis, 'edgecolor','none')
legend('background')
c = gray(256);
c = flipud(c);
colormap(c)
shading interp
colorbar
view(0,90)

%%% Test 1: assimilation every step.
y_sets = 5;
y = zeros(n, K-1, y_sets);

parfor i = 2:K+1
    for j=1:y_sets
        y(:, i-1, j) = H(i, x(:,i), R);
    end
end

min = 10;
max = 30;

[ ARMSE_b, RMSE_b ] = averageRootMeanSquareError(x_b, x);
ARMSE_MLETKF = zeros(y_sets, max);

ens = min:max;
for i=ens
    x_0b_en = ensembleInit(x_0b, P_0b, i);
    parfor j=1:y_sets
        y_cur = y(:,:,j);
        [ x_a_MLETKF, x_a_MLETKF_en] = da_hyb_MaximumLikelihoodEnsembleTransformKalmanFilter(x_0b_en, y_cur, @(t,x) M(t,x,F,T), @(t,x) H(t, x, 0*eye(n)), R, 1:K+1, 'MaxIter', 20);
        [ x_a_MLETKFInfl, x_a_MLETKFInfl_en] = da_hyb_MaximumLikelihoodEnsembleTransformKalmanFilter(x_0b_en, y_cur, @(t,x) M(t,x,F,T), @(t,x) H(t, x, 0*eye(n)), R, 1:K+1, 'Inflation', 1.03, 'MaxIter', 20);
        
        [ ARMSE_MLETKF(j, i), RMSE_MLETKF] = averageRootMeanSquareError(x_a_MLETKF, x);
        [ ARMSE_MLETKFInfl(j, i), RMSE_MLETKFInfl] = averageRootMeanSquareError(x_a_MLETKFInfl, x);
    end
    i
    save Lorenz95Test
end

ARMSE_MLETKF = mean(ARMSE_MLETKF);
ARMSE_MLETKFInfl = mean(ARMSE_MLETKFInfl);

figure
plot(0, ARMSE_b, 'k*', ens, ARMSE_MLETKF(ens), 'g*--', ens, ARMSE_MLETKFInfl(ens), 'g*-');
legend('background ARMSE', 'MLETKF ARMSE', 'MLETKF ARMSE Inflation');

time = toc

cd ../test
save Lorenz95Test