clear;
close all;

cd ../dev/

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
J_H = eye(n);

%%% RUN the model, to get the true result (given the discretization).
%   Use the T = 0.05;

r = 5;
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
y_sets = 4;
y = zeros(n, K-1, y_sets);

parfor i = 2:K+1
    for j=1:y_sets
        y(:, i-1, j) = H(i, x(:,i), R);
    end
end

min = 10;
max = 30;

c_min = 2;
c_max = 6;

[ ARMSE_b, RMSE_b ] = averageRootMeanSquareError(x_b, x);
ARMSE_ETKF = zeros(y_sets, max);
ARMSE_ETKFLoc = zeros(y_sets, max, c_max);
ARMSE_ETKFInflLoc = zeros(y_sets, max, c_max);



ens = min:max;
for i=ens
    x_0b_en = ensembleInit(x_0b, P_0b, i);
    parfor j=1:y_sets
        y_cur = y(:,:,j);
        [ x_a_ETKF, x_a_ETKF_en] = da_seq_EnsembleTransformKalmanFilter(x_0b_en, y_cur, @(t,x) M(t,x,F,T), @(t,x) H(t, x, 0*eye(n)), R);
        [ ARMSE_ETKF(j, i), RMSE_ETKF] = averageRootMeanSquareError(x_a_ETKF, x);
        
        for c= c_min:c_max
            [ x_a_ETKFLoc, x_a_ETKFLoc_en] = da_seq_LocalizedEnsembleTransformKalmanFilter(x_0b_en, y_cur, @(t,x) M(t,x,F,T), J_H, R, rho(:,:,c));
            [ x_a_ETKFInflLoc, x_a_ETKFInflLoc_en] = da_seq_LocalizedEnsembleTransformKalmanFilter(x_0b_en, y_cur, @(t,x) M(t,x,F,T), J_H, R, rho(:,:,c), 1:K+1, 'Inflation', 1.025);
            
            [ ARMSE_ETKFLoc(j, i, c), RMSE_ETKFLoc] = averageRootMeanSquareError(x_a_ETKFLoc, x);
            [ ARMSE_ETKFInflLoc(j, i, c), RMSE_ETKFInflLoc] = averageRootMeanSquareError(x_a_ETKFInflLoc, x);
        end
    end
    i
    save Lorenz95Test
end

ARMSE_ETKF = mean(ARMSE_ETKF);
ARMSE_ETKFLoc = squeeze(mean(ARMSE_ETKFLoc));
ARMSE_ETKFInflLoc = squeeze(mean(ARMSE_ETKFInflLoc));

figure
plot(0, ARMSE_b, 'k*', ens, ARMSE_ETKF(ens), 'g*--', ens, ARMSE_ETKFLoc(ens, c_min), 'c*-', ens, ARMSE_ETKFLoc(ens, c_max), 'c*--', ens, ARMSE_ETKFInflLoc(ens, c_min), 'b*-', ens, ARMSE_ETKFInflLoc(ens, c_max), 'b*--');
legend('background ARMSE', 'ETKF ARMSE', ['ETKF ARMSE Localization ' num2str(c_min)], ['ETKF ARMSE Localization ' num2str(c_max)], ['ETKF ARMSE Inflation Localization ' num2str(c_min)], ['ETKF ARMSE Inflation Localization ' num2str(c_max)]);

figure
c = c_min:c_max;
plot(c, ARMSE_ETKFLoc(max,c), c, ARMSE_ETKFInflLoc(max,c));
legend(['ARMSE ensemble ' num2str(max)], ['ARMSE ensemble ' num2str(max) 'Localization and Inflation']);

time = toc

cd ../test
save Lorenz95Test