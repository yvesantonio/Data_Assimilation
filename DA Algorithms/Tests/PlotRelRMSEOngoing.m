% Plot Comparisons
clear
close all

load ETKFTest.mat
ETKF_burn_inRelRMSE = cell(1, max);
ETKFInfl_burn_inRelRMSE = cell(1, max);
ETKFLoc_burn_inRelRMSE = cell(1, max);
ETKFLocInfl_burn_inRelRMSE = cell(1, max);

for i=ens
    [ a, ETKF_burn_inRelRMSE{i} ] = averageRelativeRootMeanSquareError(x_a{i}, x);
end

load ETKFInflTest.mat
for i=ens
    [ a, ETKFInfl_burn_inRelRMSE{i} ] = averageRelativeRootMeanSquareError(x_a{i}, x);
end

load ETKFLocTest.mat
for i=ens
    [ a, ETKFLoc_burn_inRelRMSE{i} ] = averageRelativeRootMeanSquareError(x_a{i}, x);
end

load ETKFLocInflTest.mat
for i=ens
    [ a, ETKFLocInfl_burn_inRelRMSE{i} ] = averageRelativeRootMeanSquareError(x_a{i}, x);
end

time = 1:500;
ens = 10
figure
plot(time, ETKF_burn_inRelRMSE{ens}(time), time, ETKFInfl_burn_inRelRMSE{ens}(time), time, ETKFLoc_burn_inRelRMSE{ens}(time), time, ETKFLocInfl_burn_inRelRMSE{ens}(time));
legend('ETKF', 'ETKF Inflation', 'ETKF Localization', 'ETKF Localization - Inflation')
title('Relative RMSE evolution - Ensemble 10')

ens = 25
figure
plot(time, ETKF_burn_inRelRMSE{ens}(time), time, ETKFInfl_burn_inRelRMSE{ens}(time), time, ETKFLoc_burn_inRelRMSE{ens}(time), time, ETKFLocInfl_burn_inRelRMSE{ens}(time));
legend('ETKF', 'ETKF Inflation', 'ETKF Localization', 'ETKF Localization - Inflation')
title('Relative RMSE evolution - Ensemble 25')

ens = 50
figure
plot(time, ETKF_burn_inRelRMSE{ens}(time), time, ETKFInfl_burn_inRelRMSE{ens}(time), time, ETKFLoc_burn_inRelRMSE{ens}(time), time, ETKFLocInfl_burn_inRelRMSE{ens}(time));
legend('ETKF', 'ETKF Inflation', 'ETKF Localization', 'ETKF Localization - Inflation')
title('Relative RMSE evolution - Ensemble 50')