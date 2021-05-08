clear;
close all;

load TestData

L_min = 1;
L_max = 50;

x_a = cell(1, L_max);
d4DVar_ARMSE = zeros(1, L_max);
d4DVar_RMSE = cell(1, L_max);
d4DVar_ARelRMSE = zeros(1, L_max);
d4DVar_RelRMSE = cell(1, L_max);

tic
parfor L=L_min:L_max
    x_a{L}(:,1) = x_0b;
    
    B = P0b;
    for i=2:K-L
        x_f = M(i, x_a{L}(:, i-1));
        x_a{L}(:, i) = da_var_4DVarStrong(x_f, y(:,i-1:i+L-1), M, J_M, H, J_H, B, R, @(f, x0) PolakRibiere(f, x0, 1e-2, 2));
        v = (x_a{L}(:, i) - x_f);
        B = v*v';
        % Inflation
        % B = 1.1*v*v';
        % Convex combination
        % B = 0.275*B + 0.625*v*v';
        % if range(B) < range(R) : B=R;
    end
    
    for i=1:L
        x_f = M(i, x_a{L}(:, K-L+i-1));
        x_a{L}(:, K-L+i) = da_var_4DVarStrong(x_f, y(:,K-L+i-2:end), M, J_M, H, J_H, B, R, @(f, x0) PolakRibiere(f, x0, 1e-2, 2));
        v = (x_a{L}(:, i) - x_f);
        B = 1.1*v*v';
    end
    
    [ d4DVar_ARMSE(L), d4DVar_RMSE{L} ] = averageRootMeanSquareError(x_a{L}(:,burn_in+1:K), x(:,burn_in+1:K));
    [ d4DVar_ARelRMSE(L), d4DVar_RelRMSE{L} ] = averageRelativeRootMeanSquareError(x_a{L}(:,burn_in+1:K), x(:,burn_in+1:K));
    L
end
time = toc

save d4DVarTest