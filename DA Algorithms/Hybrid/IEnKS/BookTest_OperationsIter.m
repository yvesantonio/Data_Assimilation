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

ens = 35;

iter = [ 50, 60, 70, 80 ];
max = 4;

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
% 
% da_method = @(x_0_en) da_seq_bundleIterativeEnsembleKalmanSmoother(x_0_en, y, M, H, R, bIEnKSOptions('OptAlg', @(f, x) GaussNewton(f, x, 1e-3, 10), 'S', 1, 'L', 0, 'Localization', sel, 'Inflation', 1.01));
da_method = @(x_0_en, iter) da_seq_bundleIterativeEnsembleKalmanSmoother(x_0_en, y, M, H, R, bIEnKSOptions('OptAlg', @(f, x) QuasiNewtonBFGS(f, x, 1e-6, iter), 'S', 1, 'L', 50, 'epsilon', 0.1)); % , 'Inflation', 1.023));

tic
parfor i=1:4
    [ x_a{i}, x_a_en ] = da_method(x_0_en{ens}, iter(i));
    [ ARMSE(i), RMSE{i} ] = averageRootMeanSquareError(x_a{i}(:,burn_in+1:K), x(:,burn_in+1:K));
    [ ARelRMSE(i), RelRMSE{i} ] = averageRelativeRootMeanSquareError(x_a{i}(:,burn_in+1:K), x(:,burn_in+1:K));
    P_last{i} = cov(x_a_en(:,:,end)');
    i
end
time = toc

IEnKSIter_ARMSE = ARMSE;
clear ARMSE
IEnKSIter_RMSE = RMSE;
clear RMSE
IEnKSIter_ARelRMSE = ARelRMSE;
clear ARelRMSE
IEnKSIter_RelRMSE = RelRMSE;
clear RelRMSE

save IEnKSIterTest