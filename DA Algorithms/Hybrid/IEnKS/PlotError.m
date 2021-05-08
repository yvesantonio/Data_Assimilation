clear
close all

load MLEFTest.mat
[ ARelRMSE, RelRMSE ] = averageRelativeRootMeanSquareError(x_a{50}, x);

load MLEFInflTest.mat
[ ARelRMSE, RelRMSE_Infl ] = averageRelativeRootMeanSquareError(x_a{50}, x);

figure('Renderer', 'painters', 'Position', [500 300 560 317])
plot(1:K, RelRMSE, 1:K, RelRMSE_Infl)
legend('Relative Erroe - MLEF', 'Relative Erroe - MLEF Inflated')
title('MLEF Error, ensemble size=50');
xlabel('Time Steps')
ylabel('Relative RMSE')