% Plot Comparisons
clear
close all

load d4DVarTestL1-50Iter15B_Simple.mat
load ETKSWINDOWTest30.mat
load IEnKSLWindowTest.mat
load IEnKSLWindowInflationTest.mat
load IEnKSMDALWindowTest.mat
load IEnKSMDALWindowInflationTest.mat

figure('Renderer', 'painters', 'Position', [500 300 560 317])
plot(1:25, d4DVar_ARelRMSE(1:25), 1:25, ETKS30_ARelRMSE(1:25), 1:25, ETKS30Infl_ARelRMSE(1:25), L(1:5), IEnKS_ARelRMSE(1:5), L(1:5), IEnKSInfl_ARelRMSE(1:5), L(1:5), IEnKS_MDA_ARelRMSE(1:5), L(1:5), IEnKSInfl_MDA_ARelRMSE(1:5));
legend('4D-Var', 'ETKS l=30', 'ETKS Inflation l=30', 'IEnKS SDA', 'IEnKS SDA Inflation', 'IEnKS MDA', 'IEnKS MDA Inflation');
title('Relative Error vs Window Length');
xlabel('Window Length L')
ylabel('Relative RMSE')

figure('Renderer', 'painters', 'Position', [500 300 560 317])
plot(1:25, d4DVar_ARelRMSE(1:25), 1:25, ETKS30_ARelRMSE(1:25), L(1:5), IEnKS_ARelRMSE(1:5), L(1:5), IEnKS_MDA_ARelRMSE(1:5));
legend('4D-Var', 'ETKS l=30', 'IEnKS SDA', 'IEnKS MDA');
title('Relative Error vs Window Length');
xlabel('Window Length L')
ylabel('Relative RMSE')

figure('Renderer', 'painters', 'Position', [500 300 560 317])
plot(1:50, d4DVar_ARelRMSE(1:50), 1:50, ETKS30Infl_ARelRMSE(1:50), L, IEnKSInfl_MDA_ARelRMSE(1:8));
legend('4D-Var', 'ETKS Inflation l=30', 'IEnKS MDA Inflation l=25');
title('Relative Error vs Window Length');
xlabel('Window Length L')
ylabel('Relative RMSE')