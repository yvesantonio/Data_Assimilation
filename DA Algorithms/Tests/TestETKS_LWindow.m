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

L_min = 1;
L_max = 50;

ARMSE = zeros(1, L_max);
RMSE = cell(1, L_max);
ARelRMSE = zeros(1, L_max);
RelRMSE = cell(1, L_max);

ARMSE_Infl = zeros(1, L_max);
RMSE_Infl = cell(1, L_max);
ARelRMSE_Infl = zeros(1, L_max);
RelRMSE_Infl = cell(1, L_max);

ens = 30;

tic

x_a_filter = da_seq_EnsembleTransformKalmanFilter(x_0_en{ens}, y, M, H, R);
[ x_a, x_a_en, x_a_L ] = da_seq_EnsembleTransformKalmanSmootherStoreALL(x_0_en{ens}, y, M, H, R, L_max);
[ x_a_Infl, x_a_en, x_a_Infl_L ] = da_seq_EnsembleTransformKalmanSmootherStoreALL(x_0_en{ens}, y, M, H, R, L_max, 1:K, 'Inflation', 1.025);

time = toc

parfor L_window=L_min:L_max
    [ ARMSE(L_window), RMSE{L_window} ] = averageRootMeanSquareError(x_a_L{L_window}(:,burn_in+1:K), x(:,burn_in+1:K));
    [ ARelRMSE(L_window), RelRMSE{L_window} ] = averageRelativeRootMeanSquareError(x_a_L{L_window}(:,burn_in+1:K), x(:,burn_in+1:K));
    
    [ ARMSE_Infl(L_window), RMSE_Infl{L_window} ] = averageRootMeanSquareError(x_a_Infl_L{L_window}(:,burn_in+1:K), x(:,burn_in+1:K));
    [ ARelRMSE_Infl(L_window), RelRMSE_Infl{L_window} ] = averageRelativeRootMeanSquareError(x_a_Infl_L{L_window}(:,burn_in+1:K), x(:,burn_in+1:K));   
    L_window
end

ETKS30_ARMSE = ARMSE;
clear ARMSE
ETKS30_RMSE = RMSE;
clear RMSE
ETKS30_ARelRMSE = ARelRMSE;
clear ARelRMSE
ETKS30_RelRMSE = RelRMSE;
clear RelRMSE

ETKS30Infl_ARMSE = ARMSE_Infl;
clear ARMSE_Infl
ETKS30Infl_RMSE = RMSE_Infl;
clear RMSE_Infl
ETKS30Infl_ARelRMSE = ARelRMSE_Infl;
clear ARelRMSE_Infl
ETKS30Infl_RelRMSE = RelRMSE_Infl;
clear RelRMSE_Infl

save ETKSWINDOWTest30