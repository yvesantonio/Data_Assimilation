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
ens = min:max;

selector = cell(1,n);

for i1=1:n
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
s = 10;
t_obs = 1:s:K;
y = y(:,t_obs);
da_method = @(x_0_en) da_seq_EnsembleTransformKalmanFilter(x_0_en, y, M, H, R, t_obs, 'Localization', sel, 'Inflation', 1.01);
K = t_obs(end);
tic
for i=ens
    [ x_a{i}, x_a_en ] = da_method(x_0_en{i});
    [ ARMSE(i), RMSE{i} ] = averageRootMeanSquareError(x_a{i}(:,burn_in+1:K), x(:,burn_in+1:K));
    [ ARelRMSE(i), RelRMSE{i} ] = averageRelativeRootMeanSquareError(x_a{i}(:,burn_in+1:K), x(:,burn_in+1:K));
    P_last{i} = cov(x_a_en(:,:,end)');
    save ETKFLocInflTestSampling10
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

save ETKFLocInflTestSampling10
