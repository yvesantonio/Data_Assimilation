% Plot Comparisons
clear
close all

load IEnKFInflTest.mat
[ ARelRMSE, RelRMSE ] = averageRelativeRootMeanSquareError(x_a{50}, x);

sampl = 1:50;
figure('Renderer', 'painters', 'Position', [500 300 560 317])
plot(sampl, RelRMSE(sampl))
legend('IEnKF Inflation=1.023')
title('Relative Error');
xlabel('time step')
ylabel('Relative RMSE')

load MLEFInflTest.mat
[ ARelRMSE, RelRMSE ] = averageRelativeRootMeanSquareError(x_a{50}, x);

sampl = 1:50;
figure('Renderer', 'painters', 'Position', [500 300 560 317])
plot(sampl, RelRMSE(sampl))
legend('MLEF Inflation=1.023')
title('Relative Error');
xlabel('time step')
ylabel('Relative RMSE')