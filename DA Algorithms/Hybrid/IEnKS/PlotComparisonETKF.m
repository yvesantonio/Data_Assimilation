% Plot Comparisons
clear
close all

load MLEFInflTest.mat
load ETKFInflTest.mat

figure('Renderer', 'painters', 'Position', [500 300 560 317])
plot(ens, MLEFInfl_ARelRMSE(ens), ens, ETKFInfl_ARelRMSE(ens));
legend('MLEF Inflation', 'ETKF Inflation');
title('MLEF Vs. ETKF');
xlabel('Ensemble size')
ylabel('Relative RMSE')