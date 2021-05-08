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

ens = 20;

L_min = 1;
L_max = 50;

ARMSE = zeros(1, L_max);
RMSE = cell(1, L_max);
ARelRMSE = zeros(1,  L_max);
RelRMSE = cell(1, L_max);

ARMSE_Infl = zeros(1, L_max);
RMSE_Infl = cell(1, L_max);
ARelRMSE_Infl = zeros(1, L_max);
RelRMSE_Infl = cell(1, L_max);

tic
[ x_a_mean, x_a_en, x_a ] = da_seq_EnsembleTransformKalmanSmootherStoreALL(x_0_en{ens}, y, M, H, R, L_max);
[ x_a_mean, x_a_en, x_a_Infl ] = da_seq_EnsembleTransformKalmanSmootherStoreALL(x_0_en{ens}, y, M, H, R, L_max, 1:K, 'Inflation', 1.05);

time = toc

parfor L_window=L_min:L_max
    [ ARMSE(L_window), RMSE{L_window} ] = averageRootMeanSquareError(x_a{L_window}(:,burn_in+1:K), x(:,burn_in+1:K));
    [ ARelRMSE(L_window), RelRMSE{L_window} ] = averageRelativeRootMeanSquareError(x_a{L_window}(:,burn_in+1:K), x(:,burn_in+1:K));
    
    [ ARMSE_Infl(L_window), RMSE_Infl{L_window} ] = averageRootMeanSquareError(x_a_Infl{L_window}(:,burn_in+1:K), x(:,burn_in+1:K));
    [ ARelRMSE_Infl(L_window), RelRMSE_Infl{L_window} ] = averageRelativeRootMeanSquareError(x_a_Infl{L_window}(:,burn_in+1:K), x(:,burn_in+1:K));
end

ETKS20_ARMSE = ARMSE;
clear ARMSE
ETKS20_RMSE = RMSE;
clear RMSE
ETKS20_ARelRMSE = ARelRMSE;
clear ARelRMSE
ETKS20_RelRMSE = RelRMSE;
clear RelRMSE

ETKS20Infl_ARMSE = ARMSE_Infl;
clear ARMSE_Infl
ETKS20Infl_RMSE = RMSE_Infl;
clear RMSE_Infl
ETKS20Infl_ARelRMSE = ARelRMSE_Infl;
clear ARelRMSE_Infl
ETKS20Infl_RelRMSE = RelRMSE_Infl;
clear RelRMSE_Infl

save ETKSWINDOWTest20-2