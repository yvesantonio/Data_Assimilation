% Plot Comparisons
clear
close all

load ETKFTest.mat
[ ARelRMSE, RelRMSE ] = averageRelativeRootMeanSquareError(x_a{50}, x);

sampl = 1:50;
figure('Renderer', 'painters', 'Position', [500 300 560 317])
plot(sampl, RelRMSE(sampl))
legend('ETKF')
title('Relative Error');
xlabel('time step')
ylabel('Relative RMSE')

% % % % % % % % % % % % % % % load ETKFInflTest.mat
% % % % % % % % % % % % % % % [ ARelRMSE, RelRMSEInfl ] = averageRelativeRootMeanSquareError(x_a{50}, x);
% % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % load ETKFLocTest.mat
% % % % % % % % % % % % % % % [ ARelRMSE, RelRMSELoc ] = averageRelativeRootMeanSquareError(x_a{50}, x);
% % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % load ETKFLocInflTest.mat
% % % % % % % % % % % % % % % [ ARelRMSE, RelRMSELocInfl ] = averageRelativeRootMeanSquareError(x_a{50}, x);

% % % % % % % % % % % % % % % load ETKSTest.mat
% % % % % % % % % % % % % % % [ ARelRMSE, RelRMSE_ETKS ] = averageRelativeRootMeanSquareError(x_a{50}, x);
% % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % figure('Renderer', 'painters', 'Position', [500 300 560 317])
% % % % % % % % % % % % % % % plot(sampl, RelRMSE_ETKS(sampl))
% % % % % % % % % % % % % % % legend('ETKS')
% % % % % % % % % % % % % % % title('Relative Error');
% % % % % % % % % % % % % % % xlabel('time step')
% % % % % % % % % % % % % % % ylabel('Relative RMSE')

load d4DVarTestL1-50Iter15B_Simple.mat
[ ARelRMSE, RelRMSE_4D ] = averageRelativeRootMeanSquareError(x_a{50}, x);

sampl = 1:500;
figure('Renderer', 'painters', 'Position', [500 300 560 317])
plot(sampl, RelRMSE_4D(sampl))
legend('4D-Var L=50')
title('Relative Error');
xlabel('time step')
ylabel('Relative RMSE')