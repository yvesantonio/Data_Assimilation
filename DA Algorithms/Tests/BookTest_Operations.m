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

min = 16;
max = 25;
ens = min:max;

selector = cell(1,n);

parfor i1=1:n
    selector{i1} = roundSelector(i1, n, 4);
end
sel = @(i) selector{i};

ARMSE = zeros(1, max);
RMSE = cell(1, max);
ARelRMSE = zeros(1, max);
RelRMSE = cell(1, max);
P_last = cell(1, max);

x_a = cell(1, max);
% x_a_en = cell(1, n);

da_method = @(x_0_en) da_seq_EnsembleTransformKalmanFilter(x_0_en, y, M, H, R, 1:K, 'Localization', sel, 'Inflation', 1.01);

tic
parfor i=ens
    [ x_a{i}, x_a_en ] = da_method(x_0_en{i});
    [ ARMSE(i), RMSE{i} ] = averageRootMeanSquareError(x_a{i}(:,burn_in+1:K), x(:,burn_in+1:K));
    [ ARelRMSE(i), RelRMSE{i} ] = averageRelativeRootMeanSquareError(x_a{i}(:,burn_in+1:K), x(:,burn_in+1:K));
    P_last{i} = cov(x_a_en(:,:,end)');
    i
end
time = toc

ETKFLocInfl_ARMSE = ARMSE;
clear ARMSE
ETKFLocInfl_RMSE = RMSE;
clear RMSE
ETKFLocInfl_ARelRMSE = ARelRMSE;
clear ARelRMSE
ETKFLocInfl_RelRMSE = RelRMSE;
clear RelRMSE

save ETKFLocInflTest2