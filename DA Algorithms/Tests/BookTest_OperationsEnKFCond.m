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

ACond2YY = zeros(1, max);
Cond2YY = cell(1, max);

ARMSE = zeros(1, max);
RMSE = cell(1, max);
ARelRMSE = zeros(1, max);
RelRMSE = cell(1, max);
P_last = cell(1, max);

x_a = cell(1, max);
% x_a_en = cell(1, n);

da_method = @(x_0_en) da_seq_StochasticEnKFCOND(x_0_en, y, M, H, R);

tic
parfor i=ens
    [ x_a{i}, x_a_en, Cond2YY{i} ] = da_method(x_0_en{i});
    ACond2YY(i) = mean(Cond2YY{i});
    [ ARMSE(i), RMSE{i} ] = averageRootMeanSquareError(x_a{i}(:,burn_in+1:K), x(:,burn_in+1:K));
    [ ARelRMSE(i), RelRMSE{i} ] = averageRelativeRootMeanSquareError(x_a{i}(:,burn_in+1:K), x(:,burn_in+1:K));
    P_last{i} = cov(x_a_en(:,:,end)');
    i
end
time = toc

save EnKFCond2