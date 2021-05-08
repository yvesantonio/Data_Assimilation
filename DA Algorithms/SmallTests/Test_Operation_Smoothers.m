clear;
close all;

load TestData

min = 2;
max = 25;
ens = min:max;

L_min = 1;
L_max = 50;

ARMSE = zeros(L_max, max);
RMSE = cell(L_max, max);
ARelRMSE = zeros(L_max, max);
RelRMSE = cell(L_max, max);

x_a_meanL = cell(1, max);

tic
parfor i=ens
    [ x_a_mean, x_a, x_a_meanL{i} ] = da_seq_EnsembleTransformKalmanSmootherStoreALL(x_0_en{i}, y, M, H, R, L_max, 1:K, 'Inflation', 1.0002);
    
    for L=L_min:L_max
        [ ARMSE(L, i), RMSE{L, i} ] = averageRootMeanSquareError(x_a_meanL{i}{L}(:,burn_in+1:K), x(:,burn_in+1:K));
        [ ARelRMSE(L, i), RelRMSE{L, i} ] = averageRelativeRootMeanSquareError(x_a_meanL{i}{L}(:,burn_in+1:K), x(:,burn_in+1:K));
    end
    i
end
time = toc

ETKSInfl_ARMSE = ARMSE;
clear ARMSE
ETKSInfl_RMSE = RMSE;
clear RMSE
ETKSInfl_ARelRMSE = ARelRMSE;
clear ARelRMSE
ETKSInfl_RelRMSE = RelRMSE;
clear RelRMSE

save ETKSInfl

% % % % % % % % % % % % % % ETKS_ARMSE = ARMSE;
% % % % % % % % % % % % % % clear ARMSE
% % % % % % % % % % % % % % ETKS_RMSE = RMSE;
% % % % % % % % % % % % % % clear RMSE
% % % % % % % % % % % % % % ETKS_ARelRMSE = ARelRMSE;
% % % % % % % % % % % % % % clear ARelRMSE
% % % % % % % % % % % % % % ETKS_RelRMSE = RelRMSE;
% % % % % % % % % % % % % % clear RelRMSE
% % % % % % % % % % % % % % 
% % % % % % % % % % % % % % save ETKS