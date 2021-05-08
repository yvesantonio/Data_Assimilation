clear;
close all;

load TestData

%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Run 1:
min = 5;
max = 20;
infl = 1:0.0001:1.02;
%%%%%%%%%%%%%% FILE GENERATED: ETKFOptInflTest0
%%%%%%%%%%%%%%

r = length(infl);

ARMSE = zeros(1, r);
RMSE = cell(1, r);
ARelRMSE = zeros(1, r);
RelRMSE = cell(1, r);

x_a = cell(1, r);

t_inst = 1:K;

da_method = @(x_0_en, inflation) da_seq_EnsembleTransformKalmanFilter(x_0_en, y, M, H, R, t_inst, 'Inflation', inflation);

tic
parfor h = 1:r    
    [ x_a{h}, x_a_en ] = da_method(x_0_en{3}, infl(h));
    [ ARMSE(h), RMSE{h} ] = averageRootMeanSquareError(x_a{h}(:,burn_in+1:K), x(:,burn_in+1:K));
    [ ARelRMSE(h), RelRMSE{h} ] = averageRelativeRootMeanSquareError(x_a{h}(:,burn_in+1:K), x(:,burn_in+1:K));
end
time = toc

ETKFInfl_ARMSE = ARMSE;
clear ARMSE
ETKFInfl_RMSE = RMSE;
clear RMSE
ETKFInfl_ARelRMSE = ARelRMSE;
clear ARelRMSE
ETKFInfl_RelRMSE = RelRMSE;
clear RelRMSE

save ETKFOptInflTest0

plot(infl, ETKFInfl_ARelRMSE)