clear
close all

load ETKSTest.mat
L_min = 1;
for i=5:46
    [ x_a_1_46{1, i}, x_a_en ] = da_seq_EnsembleTransformKalmanSmoother(x_0_en{i}, y, M, H, R, 1);
    [ ETKS_ARMSE(1, i), ETKS_RMSE{1, i} ] = averageRootMeanSquareError(x_a_1_46{1, i}(:,burn_in+1:K), x(:,burn_in+1:K));
    [ ETKS_ARelRMSE(1, i), ETKS_RelRMSE{1, i} ] = averageRelativeRootMeanSquareError(x_a_1_46{1, i}(:,burn_in+1:K), x(:,burn_in+1:K));
    P_last{1, i} = cov(x_a_en(:,:,end)');
    i
end
time = toc

save ETKSTest1