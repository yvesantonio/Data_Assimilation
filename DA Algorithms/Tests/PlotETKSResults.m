% Plot Comparisons
clear
close all

load ETKSTest.mat
load ETKSInflTest.mat

legendArg = cell(1, L_max);

figure('Renderer', 'painters', 'Position', [500 300 560 317])
for L=L_min:L_max
    hold on
    plot(ens, ETKS_ARelRMSE(L, ens));
    legendArg{L} = ['ETKS (L ' num2str(L) ')' ];
end

legend(legendArg);
title('Relative Error - ETKS');

figure('Renderer', 'painters', 'Position', [500 300 560 317])
for L=L_min:L_max
    hold on
    plot(ens, ETKSInfl_ARelRMSE(L, ens));
    legendArg{L} = ['ETKS Inflation (L ' num2str(L) ')' ];
end

legend(legendArg);
title('Relative Error - ETKS Inflation');

figure('Renderer', 'painters', 'Position', [500 300 560 317])
plot(ens, ETKS_ARelRMSE(1, ens), ens, ETKSInfl_ARelRMSE(1, ens), ens, ETKS_ARelRMSE(10, ens), ens, ETKSInfl_ARelRMSE(10, ens));
legend('ETKS, L=1', 'ETKS, L=1, Inflation', 'ETKS, L=10', 'ETKS, L=10, Inflation')
title('Relative Error - ETKS');
xlabel('Ensemble size')
ylabel('Relative RMSE')
