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
L = [1, 5, 10, 15, 25, 30, 40, 50];
selector = cell(1,n);

for i1=1:n
    selector{i1} = roundSelector(i1, n, 4);
end
sel = @(i) selector{i};

ARMSE = zeros(1, 50);
RMSE = cell(1, 50);
ARelRMSE = zeros(1, 50);
RelRMSE = cell(1, 50);
P_last = cell(1, 50);

x_a = cell(1, 50);
% x_a_en = cell(1, n);
% 
da_method = @(x_0_en, L) da_seq_bundleIterativeEnsembleKalmanSmoother(x_0_en, y, M, H, R, bIEnKSOptions('OptAlg', @(f, x) QuasiNewtonBFGS(f, x, 1e-9, 50), 'S', 1, 'L', L, 'epsilon', 0.01));%, 'Inflation', 1.023));

tic
parfor i=1:length(L)
    [ x_a{i}, x_a_en ] = da_method(x_0_en{ens}, L(i));
    [ ARMSE(i), RMSE{i} ] = averageRootMeanSquareError(x_a{i}(:,burn_in+1:K), x(:,burn_in+1:K));
    [ ARelRMSE(i), RelRMSE{i} ] = averageRelativeRootMeanSquareError(x_a{i}(:,burn_in+1:K), x(:,burn_in+1:K));
    P_last{i} = cov(x_a_en(:,:,end)');
    ARelRMSE(i)
    L(i)
end
time = toc

% save IEnKSLWindowInflationTest

% % % % for i=length(L):-1:2
% % % %     ARMSE(L(i)) = ARMSE(i);
% % % %     ARMSE(i) = 0;
% % % %     RMSE{L(i)} = RMSE{i};
% % % %     RMSE{i} = [];
% % % %     
% % % %     ARelRMSE(L(i)) = ARelRMSE(i);
% % % %     ARelRMSE(i) = 0;
% % % %     RelRMSE{L(i)} = RelRMSE{i};
% % % %     RelRMSE{i} = [];
% % % % end


IEnKS_ARMSE = ARMSE;
clear ARMSE
IEnKS_RMSE = RMSE;
clear RMSE
IEnKS_ARelRMSE = ARelRMSE;
clear ARelRMSE
IEnKS_RelRMSE = RelRMSE;
clear RelRMSE

save IEnKSLWindowTest