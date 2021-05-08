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

parfor i1=1:n
    selector{i1} = roundSelector(i1, n, 4);
end
sel = @(i) selector{i};

IenKFLoc_ARMSE = zeros(1, max);
IenKFLoc_RMSE = cell(1, max);
IenKFLoc_ARelRMSE = zeros(1, max);
IenKFLoc_RelRMSE = cell(1, max);
IenKFLoc_P_last = cell(1, max);

MLEV_ARMSE = zeros(1, max);
MLEV_RMSE = cell(1, max);
MLEV_ARelRMSE = zeros(1, max);
MLEV_RelRMSE = cell(1, max);
MLEV_P_last = cell(1, max);


x_a = cell(1, max);
% x_a_en = cell(1, n);
% 
MLEF = @(x_0_en) da_seq_bundleIterativeEnsembleKalmanSmoother(x_0_en, y, M, H, R, bIEnKSOptions('OptAlg', @(f, x) GaussNewton(f, x, 1e-3, 20), 'S', 1, 'L', 0, 'Localization', sel));
IEnKF = @(x_0_en) da_seq_bundleIterativeEnsembleKalmanSmoother(x_0_en, y, M, H, R, bIEnKSOptions('OptAlg', @(f, x) GaussNewton(f, x, 1e-3, 20), 'S', 1, 'L', 1, 'Localization', sel));

tic
for j=1:4
    i = 35 + (j-1)*5;
    [ x_a{i}, x_a_en ] = IEnKF(x_0_en{i});
    [ IenKFLoc_ARMSE(i), IenKFLoc_RMSE{i} ] = averageRootMeanSquareError(x_a{i}(:,burn_in+1:K), x(:,burn_in+1:K));
    [ IenKFLoc_ARelRMSE(i), IenKFLoc_RelRMSE{i} ] = averageRelativeRootMeanSquareError(x_a{i}(:,burn_in+1:K), x(:,burn_in+1:K));
    IenKFLoc_P_last{i} = cov(x_a_en(:,:,end)');
    
    [ x_a{i}, x_a_en ] = MLEV(x_0_en{i});
    [ MLEV_ARMSE(i), MLEV_RMSE{i} ] = averageRootMeanSquareError(x_a{i}(:,burn_in+1:K), x(:,burn_in+1:K));
    [ MLEV_ARelRMSE(i), MLEV_RelRMSE{i} ] = averageRelativeRootMeanSquareError(x_a{i}(:,burn_in+1:K), x(:,burn_in+1:K));
    MLEV_P_last{i} = cov(x_a_en(:,:,end)');
    i
end
time = toc


save LocTest20Iter