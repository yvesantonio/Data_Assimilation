% Plot Comparisons
clear
close all

load d4DVarTestL1-50Iter10B_Simple
d4DVar_ARelRMSE_Simple = d4DVar_ARelRMSE;

load d4DVarTestL1-50Iter10B_Inflated
d4DVar_ARelRMSE_Inflated = d4DVar_ARelRMSE;

load d4DVarTestL1-50Iter10B_Limited
d4DVar_ARelRMSE_Limited = d4DVar_ARelRMSE;

L = 1:50;
figure('Renderer', 'painters', 'Position', [500 300 560 317])
plot(L, d4DVar_ARelRMSE_Simple, L, d4DVar_ARelRMSE_Inflated , L, d4DVar_ARelRMSE_Limited);
legend('4D-Var', '4D-Var B Inflated', '4D-Var B Limited');
title('Relative Error vs Window Length (4D-Var)');
xlabel('Window L');
ylabel('Relative RMSE');

load d4DVarTestL1-50Iter15B_Simple
d4DVar_ARelRMSE_Simple = d4DVar_ARelRMSE;

load d4DVarTestL1-50Iter15B_Inflated
d4DVar_ARelRMSE_Inflated = d4DVar_ARelRMSE;

load d4DVarTestL1-50Iter15B_Limited
d4DVar_ARelRMSE_Limited = d4DVar_ARelRMSE;

L = 1:50;
figure('Renderer', 'painters', 'Position', [500 300 560 317])
plot(L, d4DVar_ARelRMSE_Simple, L, d4DVar_ARelRMSE_Inflated , L, d4DVar_ARelRMSE_Limited);
legend('4D-Var', '4D-Var B Inflated', '4D-Var B Limited');
title('Relative Error vs Window Length (4D-Var)');
xlabel('Window L');
ylabel('Relative RMSE');

load d4DVarTestL50Iter.mat
S = 5:40;
figure('Renderer', 'painters', 'Position', [500 300 560 317])
plot(S, d4DVar_ARelRMSE(S), '*-');
legend('4D-Var L=50');
title('Relative Error vs Available Measurements (4D-Var)');
xlabel('Max Iterations');
ylabel('Relative RMSE');

load d4DVarTestSamplings.mat
S = sampl;
figure('Renderer', 'painters', 'Position', [500 300 560 317])
plot([1 S], [d4DVar_ARelRMSE(40) d4DVar_ARelRMSE_Samplings(S)], '*-');
legend('4D-Var L=50, Shift S');
title('Relative Error vs Available Measurements (4D-Var)');
xlabel('Shift S');
ylabel('Relative RMSE');
