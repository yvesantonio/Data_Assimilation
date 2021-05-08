clear
close all

load EnKF
load EnKFInfl

figure
plot(ens, EnKF_ARelRMSE(ens), ens, EnKFInfl_ARelRMSE(ens))
legend('EnKF', 'EnKF + Inflation')


load ETKF
load ETKFInfl
load ETKS
load ETKSInfl

ens = 3:9;
figure
plot(ens, ETKF_ARelRMSE(ens), ens, ETKFInfl_ARelRMSE(ens), ens, ETKS_ARelRMSE(50, ens), ens, ETKSInfl_ARelRMSE(50, ens))
legend('ETKF', 'ETKF + Inflation', 'ETKS L=50', 'ETKS + Inflation L=50')

ens = 4;
l = 1:50;
figure
plot(l, ETKS_ARelRMSE(:, ens), l, ETKSInfl_ARelRMSE(:, ens))
legend('ETKS Ensemble 4', 'ETKS + Inflation Ensemble = 4')

load MLEV
load IEnKF
load IEnKSL50
load MLEVInfl
load IEnKFInfl
load IEnKSL50Infl

ens = 4:20;
figure
plot(ens, MLEV_ARelRMSE(ens), ens, IEnKF_ARelRMSE(ens), ens, IEnKS_ARelRMSE(ens))
legend('MLEV', 'IEnKF', 'IEnKS L=50')

ens = 4:20;
figure
plot(ens, MLEVInfl_ARelRMSE(ens), ens, IEnKFInfl_ARelRMSE(ens), ens, IEnKSInfl_ARelRMSE(ens))
legend('MLEV Inflation', 'IEnKF Inflation', 'IEnKS Inflation L=50')