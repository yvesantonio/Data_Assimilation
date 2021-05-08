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

d3DVar = @(x_f, y, B) da_var_4DVarStrong(x_f, y, M, J_M, H, J_H, B, R, @(f, x0) PolakRibiere(f, x0, 1e-3, 15));

tic
x_a(:,1) = x_0b;
    
B = P_0b;
for i=2:K
    x_f = M(i, x_a(:, i-1));
    x_a(:, i) = d3DVar(x_f, y(:,i-1), B);
    % JMa = J_M(i, x_a(:, i));
    % B = (0.9*B^-1 + R^-1)^-1;
    B = (B + R)/2;
    i
end

[ d3DVar_ARMSE, d3DVar_RMSE ] = averageRootMeanSquareError(x_a(:,burn_in+1:K), x(:,burn_in+1:K));
[ d3DVar_ARelRMSE, d3DVar_RelRMSE ] = averageRelativeRootMeanSquareError(x_a(:,burn_in+1:K), x(:,burn_in+1:K));

time = toc

save d3DVarTestIter15-1