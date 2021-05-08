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

min = 5;
max = 50;

L_min = 1;
L_max = 10;

ens = min:max;

ARMSE = zeros(L_max, max);
RMSE = cell(L_max, max);
ARelRMSE = zeros(L_max, max);
RelRMSE = cell(L_max, max);

ARMSE_Infl = zeros(L_max, max);
RMSE_Infl = cell(L_max, max);
ARelRMSE_Infl = zeros(L_max, max);
RelRMSE_Infl = cell(L_max, max);

x_a = cell(1, max);
x_a_Infl = cell(1, max);
% x_a_en = cell(1, max);

tic
parfor i=ens
    [ x_a_mean, x_a_en, x_a{i} ] = da_seq_EnsembleTransformKalmanSmootherStoreALL(x_0_en{i}, y, M, H, R, L_max);
    [ x_a_mean, x_a_en, x_a_Infl{i} ] = da_seq_EnsembleTransformKalmanSmootherStoreALL(x_0_en{i}, y, M, H, R, L_max, 1:K, 'Inflation', 1.025);
    i
end
time = toc

parfor j=ens
    for L_window=L_min:L_max
        [ ARMSE(L_window, j), RMSE{L_window, j} ] = averageRootMeanSquareError(x_a{j}{L_window}(:,burn_in+1:K), x(:,burn_in+1:K));
        [ ARelRMSE(L_window, j), RelRMSE{L_window, j} ] = averageRelativeRootMeanSquareError(x_a{j}{L_window}(:,burn_in+1:K), x(:,burn_in+1:K));
        
        [ ARMSE_Infl(L_window, j), RMSE_Infl{L_window, j} ] = averageRootMeanSquareError(x_a_Infl{j}{L_window}(:,burn_in+1:K), x(:,burn_in+1:K));
        [ ARelRMSE_Infl(L_window, j), RelRMSE_Infl{L_window, j} ] = averageRelativeRootMeanSquareError(x_a_Infl{j}{L_window}(:,burn_in+1:K), x(:,burn_in+1:K));
    end
    j
end

ETKS_ARMSE = ARMSE;
clear ARMSE
ETKS_RMSE = RMSE;
clear RMSE
ETKS_ARelRMSE = ARelRMSE;
clear ARelRMSE
ETKS_RelRMSE = RelRMSE;
clear RelRMSE

ETKSInfl_ARMSE = ARMSE_Infl;
clear ARMSE_Infl
ETKSInfl_RMSE = RMSE_Infl;
clear RMSE_Infl
ETKSInfl_ARelRMSE = ARelRMSE_Infl;
clear ARelRMSE_Infl
ETKSInfl_RelRMSE = RelRMSE_Infl;
clear RelRMSE_Infl

save ETKSTest