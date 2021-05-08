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

min_loc = 5;
max_loc = 12;

ens = min:max;
loc = min_loc:max_loc;

ARMSE = zeros(max_loc, max);
RMSE = cell(max_loc, max);
ARelRMSE = zeros(max_loc, max);
RelRMSE = cell(max_loc, max);
P_last = cell(max_loc, max);

x_a = cell(max_loc, max);
% x_a_en = cell(1, n);

% da_method = @(x_0_en, rho) da_seq_LocalizedStochasticEnsembleKalmanFilter(x_0_en, y, M, eye(n), R, rho); % 
da_method = @(x_0_en, rho) da_seq_LocalizedEnsembleTransformKalmanFilter(x_0_en, y, M, eye(n), R, rho); % 
tic
for l=loc
    rho = roundRho(n, l);
    parfor i=ens
        [ x_a{l, i}, x_a_en ] = da_method(x_0_en{i}, rho);
        [ ARMSE(l, i), RMSE{l, i} ] = averageRootMeanSquareError(x_a{l, i}(:,burn_in+1:K), x(:,burn_in+1:K));
        [ ARelRMSE(l, i), RelRMSE{l, i} ] = averageRelativeRootMeanSquareError(x_a{l, i}(:,burn_in+1:K), x(:,burn_in+1:K));
        P_last{l, i} = cov(x_a_en(:,:,end)');
        i
    end
    l
end
time = toc

CLocETKF_ARMSE = ARMSE;
clear ARMSE
CLocETKF_RMSE = RMSE;
clear RMSE
CLocETKF_ARelRMSE = ARelRMSE;
clear ARelRMSE
CLocETKF_RelRMSE = RelRMSE;
clear RelRMSE

save CLocETKFTest