% Plot Comparisons
clear
close all

load ETKFTest.mat
load ETKFLocInflTest.mat
load ETKSWINDOWTest25.mat
load ETKSWINDOWTest30.mat
load d4DVarTestL1-50Iter15B_Simple.mat

MEASUREMENTS_ARelRMSE = averageRelativeRootMeanSquareError(y, x(:,2:end));

L = 1:50;
figure('Renderer', 'painters', 'Position', [500 300 560 317])
plot(L, ETKF_ARelRMSE(30).*ones(size(L)), L, ETKFLocInfl_ARelRMSE(10).*ones(size(L)), L, d4DVar_ARelRMSE, L, ETKS30_ARelRMSE, L, ETKS25Infl_ARelRMSE);
legend('ETKF, l=30', 'ETKF Localization Inflation, l=10', '4D-Var', 'ETKS, l=30', 'ETKS Inflation, l=25');
title(['Relative Error vs Window Length, Measurements Accuracy: ' num2str(MEASUREMENTS_ARelRMSE*100, '%1.0f') '%' ]);
xlabel('Window length')
ylabel('Relative RMSE')
xlim([1 50])
ylim([d4DVar_ARelRMSE(50)-0.01 d4DVar_ARelRMSE(1)+0.01])

% % % figure('Renderer', 'painters', 'Position', [500 300 560 317])
% % % plot(ens, ETKF_ARMSE(ens), ens, ETKFLocInfl_ARMSE(ens), ens, ETKS_ARMSE(1, ens), ens, ETKS_ARMSE(10, ens), ens, ETKSInfl_ARelRMSE(1, ens), ens, ETKSInfl_ARelRMSE(10, ens));
% % % legend('ETKF', 'ETKF Localization Inflation', 'ETKS L=2', 'ETKS L=10', 'ETKS Inflation L=1', 'ETKS Inflation L=10');
% % % title('ETKF vs ETKS');
% % % xlabel('Ensemble size')
% % % ylabel('Relative RMSE')
