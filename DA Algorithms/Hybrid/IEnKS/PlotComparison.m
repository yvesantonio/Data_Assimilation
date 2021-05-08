% Plot Comparisons
clear
close all

load MLEFTest.mat
load MLEFInflTest.mat
load MLEFLocTest.mat
load MLEFLocInflTest.mat

figure('Renderer', 'painters', 'Position', [500 300 560 317])
plot(ens, MLEF_ARelRMSE(ens), ens, MLEFInfl_ARelRMSE(ens), ens, MLEFLoc_ARelRMSE(ens), ens, MLEFLocInfl_ARelRMSE(ens));
legend('MLEF', 'MLEF Inflation', 'MLEF Localization', 'MLEF Localization Inflation');
title('Relative Error vs Ensemble Size (MLEF)');
xlabel('Ensemble size')
ylabel('Relative RMSE')

clear

load IEnKFTest.mat
load IEnKFInflTest.mat
load IEnKFLocTest.mat
load IEnKFLocInflTest.mat

figure('Renderer', 'painters', 'Position', [500 300 560 317])
plot(ens, IEnKF_ARelRMSE(ens), ens, IEnKFInfl_ARelRMSE(ens), ens, IEnKFLoc_ARelRMSE(ens), ens, IEnKFLocInfl_ARelRMSE(ens));
legend('IEnKF', 'IEnKF Inflation', 'IEnKF Localization', 'IEnKF Localization Inflation');
title('Relative Error vs Ensemble Size (IEnKF Smoothing performance)');
xlabel('Ensemble size')
ylabel('Relative RMSE')



% % % figure('Renderer', 'painters', 'Position', [500 300 560 317])
% % % plot(ens, ETKF_ARMSE(ens), ens, ETKFLocInfl_ARMSE(ens), ens, ETKS_ARMSE(1, ens), ens, ETKS_ARMSE(10, ens), ens, ETKSInfl_ARelRMSE(1, ens), ens, ETKSInfl_ARelRMSE(10, ens));
% % % legend('ETKF', 'ETKF Localization Inflation', 'ETKS L=2', 'ETKS L=10', 'ETKS Inflation L=1', 'ETKS Inflation L=10');
% % % title('ETKF vs ETKS');
% % % xlabel('Ensemble size')
% % % ylabel('Relative RMSE')
