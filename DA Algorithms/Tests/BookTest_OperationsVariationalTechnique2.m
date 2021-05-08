%%%
%%% Book Test: OPERATIONS
%%%

%%%
%%% WARNING: Make sure to change the names of the variables and save the
%%% results of the script.
%%%

clear;
close all;

load BookTestData

I = eye(n);
J_H = @(i, x) I;

%%% Compute J_M with finite differences (simplest rule: discretize, then
%%% adjoint).

% J_M = @(i, x) lorentz95Adjoint(i, x, T);
J_M = @(i, x) finiteDifferenceAdjoint(i, x, M);

x_a = zeros(n, K);

L = 50;
flag = zeros(1, length(y(1,:)));

tic
x_a(:,1) = x_0b;
x_f = M(i, x_a(:, 1));
x_a(:, i) = da_var_4DVarStrong(x_f, y(:, 1:L-1), M, J_M, H, J_H, P_0b, R, @(f, x0) PolakRibiere(f, x0, 1e-2, 10));
v = (x_a(:, 1) - x_f);
B = v*v';

for i=3:K-L
    x_f = M(i, x_a(:, i-1));
    x_a(:, i) = da_var_4DVarStrong(x_f, y(:, i+L-1), M, J_M, H, J_H, B, R, @(f, x0) PolakRibiere(f, x0, 1e-2, 20), L);
    v = (x_a(:, i) - x_f);
    B = v*v';
    i
end

for i=1:L
    x_a(:, K-L+i) = M(i, x_a(:, K-L+i-1));
end
time = toc

[ d4DVar_ARMSE, d4DVar_RMSE ] = averageRootMeanSquareError(x_a(:,burn_in+1:K), x(:,burn_in+1:K));
[ d4DVar_ARelRMSE, d4DVar_RelRMSE ] = averageRelativeRootMeanSquareError(x_a(:,burn_in+1:K), x(:,burn_in+1:K));

save d4DVarTest2L50PR80