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

L = 50;

max_sampl = 25;
% sampl = [ 1 2 5 10 25];
sampl = [ 2 5 10 25];
x_a = cell(1, max_sampl);
d4DVar_ARMSE = zeros(1, max_sampl);
d4DVar_RMSE = cell(1, max_sampl);
d4DVar_ARelRMSE = zeros(1, max_sampl);
d4DVar_RelRMSE = cell(1, max_sampl);

srR = max(eig(R));

tic
for s=sampl
    x_a{s} = zeros(n, K);
    x_a{s}(:,1) = x_0b;
    
    t_obs = 1:s:L+1;
    
    B = P_0b;
    for i=2:s:K-L
        x_f = M(i, x_a{s}(:, i-1));
        % i-1+t_obs
        x_a{s}(:, i) = da_var_4DVarStrong(x_f, y(:,i-2+t_obs), M, J_M, H, J_H, B, R, @(f, x0) PolakRibiere(f, x0, 1e-2, 40), t_obs);
        v = (x_a{s}(:, i) - x_f);
        B = v*v';
        % Inflation
        % B = 1.1*(v*v');
        % Convex combination
        %if max(eig(B)) < (1/2)*srR
        %    B = (1/2)*R;
        %end
        
        for j=1:s-1
            x_a{s}(:, i+j) = M(i, x_a{s}(:, i+j-1));
            % B = J_M(i+j+1, x_a{s}(:, i+j))'*B*J_M(i+j+1, x_a{s}(:, i+j));
        end
    end
    
    for i=0:s:L
        y_obs = K-L+i-1:s:K-1;
        t_obs = 1:s:length(y_obs)*s;
        x_f = M(K-L+i, x_a{s}(:, K-L+i-1));
        x_a{s}(:, K-L+i-1) = da_var_4DVarStrong(x_f, y(:,y_obs), M, J_M, H, J_H, B, R, @(f, x0) PolakRibiere(f, x0, 1e-2, 40), t_obs);
        v = (x_a{s}(:, K-L+i-1) - x_f);
        B = v*v';
        % B = 1.1*(v*v');
        % B = 0.125*B + 0.875*(v*v');
        % srB(K-L+i) = max(eig(B));
        for j=1:s
            x_a{s}(:, K-L+i+j) = M(K-L+i+j, x_a{s}(:, K-L+i+j-1));
            % B = J_M(K-L+i+j+1, x_a{s}(:, K-L+i+j))'*B*J_M(K-L+i+j+1, x_a{s}(:, K-L+i+j));
        end
    end
    
    [ d4DVar_ARMSE(s), d4DVar_RMSE{s} ] = averageRootMeanSquareError(x_a{s}(:,burn_in+1:K), x(:,burn_in+1:K))
    [ d4DVar_ARelRMSE(s), d4DVar_RelRMSE{s} ] = averageRelativeRootMeanSquareError(x_a{s}(:,burn_in+1:K), x(:,burn_in+1:K))
    s
end
time = toc

d4DVar_ARMSE_Samplings = d4DVar_ARMSE;
d4DVar_RMSE_Samplings = d4DVar_RMSE;
d4DVar_ARelRMSE_Samplings = d4DVar_ARelRMSE;
d4DVar_RelRMSE_Samplings = d4DVar_RelRMSE;

save d4DVarTestSamplingNoProp