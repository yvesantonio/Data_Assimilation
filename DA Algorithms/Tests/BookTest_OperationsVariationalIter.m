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

J_M = @(i, x) finiteDifferenceAdjoint(i, x, M);

L = 50;

min_iter = 5;
max_iter = 40;

x_a = cell(1, max_iter);
d4DVar_ARMSE = zeros(1, max_iter);
d4DVar_RMSE = cell(1, max_iter);
d4DVar_ARelRMSE = zeros(1, max_iter);
d4DVar_RelRMSE = cell(1, max_iter);

d4DVar = @(x_f, y, B, iter) da_var_4DVarStrong(x_f, y, M, J_M, H, J_H, B, R, @(f, x0) PolakRibiere(f, x0, 1e-2, iter));

tic
parfor iter=min_iter:max_iter
    x_a{iter}(:,1) = x_0b;
    
    B = P_0b;
    for i=2:K-L
        x_f = M(i, x_a{iter}(:, i-1));
        x_a{iter}(:, i) = d4DVar(x_f, y(:,i-1:i+L-1), B, iter);
        v = (x_a{iter}(:, i) - x_f);
        B = v*v';
    end
    
    for i=1:L
        x_f = M(i, x_a{iter}(:, K-L+i-1));
        x_a{iter}(:, K-L+i) = d4DVar(x_f, y(:,K-L+i-2:end), B, n);
        v = (x_a{iter}(:, i) - x_f);
        B = v*v';
    end
    
    [ d4DVar_ARMSE(iter), d4DVar_RMSE{iter} ] = averageRootMeanSquareError(x_a{iter}(:,burn_in+1:K), x(:,burn_in+1:K));
    [ d4DVar_ARelRMSE(iter), d4DVar_RelRMSE{iter} ] = averageRelativeRootMeanSquareError(x_a{iter}(:,burn_in+1:K), x(:,burn_in+1:K));
    iter
    % save d4DVarTestL50Iter
end
time = toc

save d4DVarTestL50Iter
