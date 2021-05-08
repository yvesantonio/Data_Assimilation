clear;
close all;

load TestData

min = 2;
max = 25;
ens = min:max;

ARMSE = zeros(1, max);
RMSE = cell(1, max);
ARelRMSE = zeros(1, max);
RelRMSE = cell(1, max);
P_last = cell(1, max);

x_a = cell(1, max);
x_a_en = cell(1, max);

% da_method = @(x_0_en) da_seq_EnsembleTransformKalmanFilter(x_0_en, y, M, H, R, 1:K, 'Inflation', 1.0002);
da_method = @(x_0_en) da_seq_StochasticEnsembleKalmanFilter(x_0_en, y, M, H, R);%, 1:K, 'Inflation', 1.0002);

tic
parfor i=ens
    [ x_a{i}, x_a_en{i} ] = da_method(x_0_en{i});
    [ ARMSE(i), RMSE{i} ] = averageRootMeanSquareError(x_a{i}(:,burn_in+1:K), x(:,burn_in+1:K));
    [ ARelRMSE(i), RelRMSE{i} ] = averageRelativeRootMeanSquareError(x_a{i}(:,burn_in+1:K), x(:,burn_in+1:K));
    P_last{i} = cov(x_a_en{i}(:,:,end)');
    i
end
time = toc

EnKF_ARMSE = ARMSE;
clear ARMSE
EnKF_RMSE = RMSE;
clear RMSE
EnKF_ARelRMSE = ARelRMSE;
clear ARelRMSE
EnKF_RelRMSE = RelRMSE;
clear RelRMSE

save EnKF