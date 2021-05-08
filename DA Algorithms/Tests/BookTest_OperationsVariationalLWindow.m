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

L_min = 1;
L_max = 50;

x_a = cell(1, L_max);
d4DVar_ARMSE = zeros(1, L_max);
d4DVar_RMSE = cell(1, L_max);
d4DVar_ARelRMSE = zeros(1, L_max);
d4DVar_RelRMSE = cell(1, L_max);

srR = max(eig(R));

tic
parfor L=L_min:L_max
    x_a{L}(:,1) = x_0b;
    
    B = P_0b;
    for i=2:K-L
        x_f = M(i, x_a{L}(:, i-1));
        x_a{L}(:, i) = da_var_4DVarStrong(x_f, y(:,i-1:i+L-1), M, J_M, H, J_H, B, R, @(f, x0) PolakRibiere(f, x0, 1e-2, 15));
        v = (x_a{L}(:, i) - x_f);
        
        % Simple
        % B = v*v';        
        
        % Inflation
        B = 1.1*(v*v');
        
        % Limited
        % if max(eig(B)) < (1/2)*srR
        %     B = (1/2)*R;
        % end
        
        % Hessian
        % He_proj = eye(n);        
        % He = B\He_proj;
        % He = 0;
        % xf = x_a{L}(:, i);
        % for j=0:L-1
        %    J_M_cur = J_M(i+j+1, xf);
        %    He_proj = J_M_cur*He_proj;
        %    xf = M(i+j+1, xf);
        %    
        %    He = He + He_proj'*J_H(i+j+1, xf)'*(R^-1)*J_H(i+j+1, xf)*He_proj;
        % end
        % B = J_M(i, x_a{L}(:, i))' * He^-1 * J_M(i, x_a{L}(:, i));
        % eig(B)
        % pause
        % i
    end
    
    for i=1:L
        x_f = M(i, x_a{L}(:, K-L+i-1));
        x_a{L}(:, K-L+i) = da_var_4DVarStrong(x_f, y(:,K-L+i-2:end), M, J_M, H, J_H, B, R, @(f, x0) PolakRibiere(f, x0, 1e-2, 40));
        v = (x_a{L}(:, i) - x_f);
        % B = v*v';
        B = 1.1*(v*v');
        % if max(eig(B)) < (1/2)*srR
        %     B = (1/2)*R;
        % end                
    end
    
    [ d4DVar_ARMSE(L), d4DVar_RMSE{L} ] = averageRootMeanSquareError(x_a{L}(:,burn_in+1:K), x(:,burn_in+1:K));
    [ d4DVar_ARelRMSE(L), d4DVar_RelRMSE{L} ] = averageRelativeRootMeanSquareError(x_a{L}(:,burn_in+1:K), x(:,burn_in+1:K));
    L
end
time = toc

save d4DVarTestL1-50Iter15B_Inflated