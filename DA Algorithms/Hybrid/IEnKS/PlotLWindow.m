load IEnKSLWindowTest.mat
load IEnKSLWindowInflationTest.mat

L_y = 1:8;
figure('Renderer', 'painters', 'Position', [500 300 560 317])
plot(L, IEnKS_ARelRMSE(L_y), L, IEnKSInfl_ARelRMSE(L_y));
legend('IEnKS', 'IEnKS Inflation');
title('Smoothing Performance SDA');
xlabel('L')
ylabel('Relative RMSE')

load IEnKSMDALWindowTest.mat
load IEnKSMDALWindowInflationTest.mat

figure('Renderer', 'painters', 'Position', [500 300 560 317])
plot(L, IEnKS_MDA_ARelRMSE(L_y), L, IEnKSInfl_MDA_ARelRMSE(L_y));
legend('IEnKS', 'IEnKS Inflation');
title('Smoothing Performance MDA');
xlabel('L')
ylabel('Relative RMSE')