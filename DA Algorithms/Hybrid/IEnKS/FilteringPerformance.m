clear;
close all;

max = 50;

x_a_filtering = cell(1, max);
% x_a_filtering_Infl = cell(1, max);

IEnKF_Filtering_ARMSE = zeros(1, max);
IEnKF_Filtering_RMSE = cell(1, max);
IEnKF_Filtering_ARelRMSE = zeros(1, max);
IEnKF_Filtering_RelRMSE = cell(1, max);

load 'IEnKFTest.mat'
parfor esize=ens
    x_a_filtering{esize} = zeros(n, K);
    for i=2:K-1
        x_a_filtering{esize}(:, i) = M(i, x_a{esize}(:,i-1));
    end
    [ IEnKF_Filtering_ARMSE(esize), IEnKF_Filtering_RMSE{esize} ]= averageRootMeanSquareError(x_a_filtering{esize}(:, burn_in+1:end), x(:, burn_in+1:end));
    [ IEnKF_Filtering_ARelRMSE(esize), IEnKF_Filtering_RelRMSE{esize} ]= averageRelativeRootMeanSquareError(x_a_filtering{esize}(:, burn_in+1:end), x(:, burn_in+1:end));
end


x_a_10filtering = cell(1, max);
% x_a_10filtering_Infl = cell(1, max);

IEnKS10_Filtering_ARMSE = zeros(1, max);
IEnKS10_Filtering_RMSE = cell(1, max);
IEnKS10_Filtering_ARelRMSE = zeros(1, max);
IEnKS10_Filtering_RelRMSE = cell(1, max);

load 'IEnKS10Test.mat'
parfor esize=ens
    x_a_10filtering{esize} = zeros(n, K);
    for i=11:K
        x_a_tmp = M(i-9, x_a{esize}(:,i-10));
        for ii=2:10
            x_a_tmp = M(i-ii+10, x_a_tmp);
        end
        x_a_10filtering{esize}(:, i) = x_a_tmp;
    end
    [ IEnKS10_Filtering_ARMSE(esize), IEnKS10_Filtering_RMSE{esize} ]= averageRootMeanSquareError(x_a_10filtering{esize}(:, burn_in+1:end), x(:, burn_in+1:end));
    [ IEnKS10_Filtering_ARelRMSE(esize), IEnKS10_Filtering_RelRMSE{esize} ]= averageRelativeRootMeanSquareError(x_a_10filtering{esize}(:, burn_in+1:end), x(:, burn_in+1:end));
end


x_a_15filtering = cell(1, max);
% x_a_10filtering_Infl = cell(1, max);

IEnKS15_Filtering_ARMSE = zeros(1, max);
IEnKS15_Filtering_RMSE = cell(1, max);
IEnKS15_Filtering_ARelRMSE = zeros(1, max);
IEnKS15_Filtering_RelRMSE = cell(1, max);

load 'IEnKS15Test.mat'
parfor esize=ens
    x_a_15filtering{esize} = zeros(n, K);
    for i=16:K
        x_a_tmp = M(i-14, x_a{esize}(:,i-15));
        for ii=2:15
            x_a_tmp = M(i-ii+15, x_a_tmp);
        end
        x_a_15filtering{esize}(:, i) = x_a_tmp;
    end
    [ IEnKS15_Filtering_ARMSE(esize), IEnKS15_Filtering_RMSE{esize} ]= averageRootMeanSquareError(x_a_15filtering{esize}(:, burn_in+1:end), x(:, burn_in+1:end));
    [ IEnKS15_Filtering_ARelRMSE(esize), IEnKS15_Filtering_RelRMSE{esize} ]= averageRelativeRootMeanSquareError(x_a_15filtering{esize}(:, burn_in+1:end), x(:, burn_in+1:end));
end


x_a_20filtering = cell(1, max);
% x_a_10filtering_Infl = cell(1, max);

IEnKS20_Filtering_ARMSE = zeros(1, max);
IEnKS20_Filtering_RMSE = cell(1, max);
IEnKS2_Filtering_ARelRMSE = zeros(1, max);
IEnKS20_Filtering_RelRMSE = cell(1, max);

load 'IEnKS20Test.mat'
parfor esize=ens
    x_a_20filtering{esize} = zeros(n, K);
    for i=21:K
        x_a_tmp = M(i-19, x_a{esize}(:,i-20));
        for ii=2:20
            x_a_tmp = M(i-ii+20, x_a_tmp);
        end
        x_a_20filtering{esize}(:, i) = x_a_tmp;
    end
    [ IEnKS20_Filtering_ARMSE(esize), IEnKS20_Filtering_RMSE{esize} ]= averageRootMeanSquareError(x_a_20filtering{esize}(:, burn_in+1:end), x(:, burn_in+1:end));
    [ IEnKS20_Filtering_ARelRMSE(esize), IEnKS20_Filtering_RelRMSE{esize} ]= averageRelativeRootMeanSquareError(x_a_20filtering{esize}(:, burn_in+1:end), x(:, burn_in+1:end));
end


figure('Renderer', 'painters', 'Position', [500 300 560 317])
plot(ens, IEnKF_Filtering_ARelRMSE(ens), ens, IEnKF_ARelRMSE(ens));
legend('IEnKF Filtering', 'IEnKF Smoothing');
title('Smoothing Performance vs Filtering performance');
xlabel('Ensemble size')
ylabel('Relative RMSE')

figure('Renderer', 'painters', 'Position', [500 300 560 317])
plot(ens, IEnKS10_Filtering_ARelRMSE(ens), ens, IEnKS10_ARelRMSE(ens))
legend('IEnKS L=10 Filtering', 'IEnKS L=10 Smoothing');
title('Smoothing Performance vs Filtering performance L=10');
xlabel('Ensemble size')
ylabel('Relative RMSE')

figure('Renderer', 'painters', 'Position', [500 300 560 317])
plot(ens, IEnKS15_Filtering_ARelRMSE(ens), ens, IEnKS15_ARelRMSE(ens))
legend('IEnKS L=15 Filtering', 'IEnKS L=15 Smoothing');
title('Smoothing Performance vs Filtering performance L=15');
xlabel('Ensemble size')
ylabel('Relative RMSE')

figure('Renderer', 'painters', 'Position', [500 300 560 317])
plot(ens, IEnKS20_Filtering_ARelRMSE(ens), ens, IEnKS20_ARelRMSE(ens))
legend('IEnKS L=20 Filtering', 'IEnKS L=20 Smoothing');
title('Smoothing Performance vs Filtering performance L=20');
xlabel('Ensemble size')
ylabel('Relative RMSE')

figure('Renderer', 'painters', 'Position', [500 300 560 317])
plot(35:50, IEnKS15_Filtering_ARelRMSE(35:50), 35:50, IEnKS15_ARelRMSE(35:50), 35:50, IEnKS20_Filtering_ARelRMSE(35:50), 35:50, IEnKS20_ARelRMSE(35:50))
legend('IEnKS L=15 Filtering', 'IEnKS L=15 Smoothing', 'IEnKS L=20 Filtering', 'IEnKS L=20 Smoothing');
title('Smoothing Performance vs Filtering performance L=20');
xlabel('Ensemble size')
ylabel('Relative RMSE')


% % % % % load 'IEnKFLocInflTest.mat'
% % % % % parfor esize=ens
% % % % %     x_a_filtering_LocInfl{esize} = zeros(n, K);
% % % % %     for i=2:K-1
% % % % %         x_a_filtering_LocInfl{esize}(:, i) = M(i, x_a{esize}(:,i-1));
% % % % %     end
% % % % %     [ IEnKF_Filtering_LocInfl_ARMSE(esize), IEnKF_LocInflFiltering_RMSE{esize} ]= averageRootMeanSquareError(x_a_filtering_LocInfl{esize}(:, burn_in+1:end), x(:, burn_in+1:end));
% % % % %     [ IEnKF_Filtering_LocInfl_ARelRMSE(esize), IEnKF_Filtering_LocInfl_RelRMSE{esize} ]= averageRelativeRootMeanSquareError(x_a_filtering_LocInfl{esize}(:, burn_in+1:end), x(:, burn_in+1:end));
% % % % % end
% % % % % 
% % % % % figure('Renderer', 'painters', 'Position', [500 300 560 317])
% % % % % plot(ens, IEnKF_Filtering_ARelRMSE(ens), ens, IEnKF_ARelRMSE(ens), ens, IEnKF_Filtering_LocInfl_ARelRMSE(ens), ens, IEnKFLocInfl_ARelRMSE(ens))
