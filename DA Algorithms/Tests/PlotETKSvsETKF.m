% Plot Comparisons
clear
close all

load ETKFTest.mat
load ETKFInflTest.mat
load ETKSWINDOWTest20.mat

ens = 20;

figure('Renderer', 'painters', 'Position', [500 300 560 317])
L = L_min:L_max;

etkf = ETKF_ARelRMSE(ens)*ones(1, L_max + 1);
etkfinfl = ETKFInfl_ARelRMSE(ens)*ones(1, L_max + 1);

plot([0 L], etkf)
plot([0 L], etkf, [0 L], [ etkf(1) ETKS20_ARelRMSE(L) ])%, [0 L], etkfinfl, [0 L], [ etkfinfl(1) ETKS20Infl_ARelRMSE(L) ]);
legend('etkf', 'etks')%, 'etkf inflation 1.025', 'etks inflation 1.025')
title('ETKF vs ETKS - Ensemble 20')
xlabel('Window Length')
ylabel('Relative RMSE')

figure('Renderer', 'painters', 'Position', [500 300 560 317])
L = L_min:L_max;

etkf = ETKF_ARelRMSE(ens)*ones(1, L_max + 1);
etkfinfl = ETKFInfl_ARelRMSE(ens)*ones(1, L_max + 1);

plot([0 L], etkf)
plot([0 L], etkfinfl, [0 L], [ etkfinfl(1) ETKS20Infl_ARelRMSE(L) ]);
legend('ETKF Inflation 1.025', 'ETKS Inflation 1.025')
title('ETKF vs ETKS - Ensemble 20')
xlabel('Window Length')
ylabel('Relative RMSE')

load ETKSWINDOWTest25.mat

ens = 25;

figure('Renderer', 'painters', 'Position', [500 300 560 317])
L = L_min:L_max;

etkf = ETKF_ARelRMSE(ens)*ones(1, L_max + 1);
etkfinfl = ETKFInfl_ARelRMSE(ens)*ones(1, L_max + 1);

plot([0 L], etkf)
plot([0 L], etkf, [0 L], [ etkf(1) ETKS25_ARelRMSE(L) ], [0 L], etkfinfl, [0 L], [ etkfinfl(1) ETKS25Infl_ARelRMSE(L) ]);
legend('etkf', 'etks', 'etkf inflation', 'etks inflation')
title('ETKF vs ETKS - Ensemble 25')
xlabel('Window Length')
ylabel('Relative RMSE')

load ETKSWINDOWTest30.mat

ens = 30;

figure('Renderer', 'painters', 'Position', [500 300 560 317])
L = L_min:L_max;

etkf = ETKF_ARelRMSE(ens)*ones(1, L_max + 1);
etkfinfl = ETKFInfl_ARelRMSE(ens)*ones(1, L_max + 1);

plot([0 L], etkf)
plot([0 L], etkf, [0 L], [ etkf(1) ETKS30_ARelRMSE(L) ], [0 L], etkfinfl, [0 L], [ etkfinfl(1) ETKS30Infl_ARelRMSE(L) ]);
legend('etkf', 'etks', 'etkf inflation', 'etks inflation')
title('ETKF vs ETKS - Ensemble 30')
xlabel('Window Length')
ylabel('Relative RMSE')