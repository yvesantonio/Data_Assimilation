% Plot Comparisons
clear
close all

load ETKFTest.mat
load ETKFInflTest.mat
load ETKFLocTest.mat
load ETKFLocInflTest.mat
load ETKSTest.mat

figure('Renderer', 'painters', 'Position', [500 300 560 317])
plot(ens, ETKF_ARelRMSE(ens), ens, ETKFInfl_ARelRMSE(ens), ens, ETKFLoc_ARelRMSE(ens), ens, ETKFLocInfl_ARelRMSE(ens));
legend('ETKF', 'ETKF Inflation', 'ETKF Localization', 'ETKF Localization Inflation');
title('Relative Error vs Ensemble Size (ETKF)');
xlabel('Ensemble size')
ylabel('Relative RMSE')

figure('Renderer', 'painters', 'Position', [500 300 560 317])
plot(ens, ETKF_ARelRMSE(ens), ens, ETKFInfl_ARelRMSE(ens), ens, ETKS_ARelRMSE(10, ens), ens, ETKSInfl_ARelRMSE(10, ens)); %, ens, ones(size(ens))*0.05,'r');
legend('ETKF', 'ETKF Inflation', 'ETKS L=10', 'ETKS Inflation L=10');
title('Relative Error vs Ensemble Size (ETKF, ETKS)');
xlabel('Ensemble size')
ylabel('Relative RMSE')

% % % figure('Renderer', 'painters', 'Position', [500 300 560 317])
% % % plot(ens, ETKF_ARMSE(ens), ens, ETKFLocInfl_ARMSE(ens), ens, ETKS_ARMSE(1, ens), ens, ETKS_ARMSE(10, ens), ens, ETKSInfl_ARelRMSE(1, ens), ens, ETKSInfl_ARelRMSE(10, ens));
% % % legend('ETKF', 'ETKF Localization Inflation', 'ETKS L=2', 'ETKS L=10', 'ETKS Inflation L=1', 'ETKS Inflation L=10');
% % % title('ETKF vs ETKS');
% % % xlabel('Ensemble size')
% % % ylabel('Relative RMSE')
