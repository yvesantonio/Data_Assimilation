%%%
%%% Inflation: look for the optimum inflation
%%%

%%%
%%% WARNING: Make sure to change the names of the variables and save the
%%% results of the script.
%%%

clear;
close all;

load BookTestData

%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Run 1:
% % % % % % % % min = 5;
% % % % % % % % max = 20;
% % % % % % % % infl = 1.01:0.05:1.31;
%%%%%%%%%%%%%% FILE GENERATED: ETKFOptInflTest0
%%%%%%%%%%%%%% 

%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Run 2:
% % % % % % % % min = 15;
% % % % % % % % max = 25;
% % % % % % % % infl = 1.01:0.02:1.1;
%%%%%%%%%%%%%% FILE GENERATED: ETKFOptInflTest1
%%%%%%%%%%%%%% 

%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Run 3:
% % % % % % % % min = 15;
% % % % % % % % max = 20;
% % % % % % % % infl = 1.01:0.005:1.04;
%%%%%%%%%%%%%% FILE GENERATED: ETKFOptInflTest2
%%%%%%%%%%%%%% 

%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Run 4:
min = 19;
max = 25;
infl = 1.005:0.005:1.035;
%%%%%%%%%%%%%% FILE GENERATED: ETKFOptInflTest3
%%%%%%%%%%%%%% 

%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Run 5:
% % % % % % % min = 17;
% % % % % % % max = 20;
% % % % % % % infl = 1.02:0.001:1.03;
%%%%%%%%%%%%%% FILE GENERATED: ETKFOptInflTest4
%%%%%%%%%%%%%% 

%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Run 6:
% % % % % % % % % min = 17;
% % % % % % % % % max = 20;
% % % % % % % % % infl = 1.03:0.001:1.035;
%%%%%%%%%%%%%% FILE GENERATED: ETKFOptInflTest5
%%%%%%%%%%%%%% 

ens = min:max;

r = length(infl);

ARMSE = zeros(r, max);
RMSE = cell(r, max);
ARelRMSE = zeros(r, max);
RelRMSE = cell(r, max);

x_a = cell(r, n);
% x_a_en = cell(r, n);

t_inst = 1:K;

da_method = @(x_0_en, inflation) da_seq_EnsembleTransformKalmanFilter(x_0_en, y, M, H, R, t_inst, 'Inflation', inflation);

tic
for h = 1:r
    parfor i=ens
        [ x_a{h, i}, x_a_en ] = da_method(x_0_en{i}, infl(h));
        [ ARMSE(h, i), RMSE{h, i} ] = averageRootMeanSquareError(x_a{h, i}(:,burn_in+1:K), x(:,burn_in+1:K));
        [ ARelRMSE(h, i), RelRMSE{h, i} ] = averageRelativeRootMeanSquareError(x_a{h, i}(:,burn_in+1:K), x(:,burn_in+1:K));
    end
    h/r
end
time = toc

ETKFInfl_ARMSE = ARMSE;
clear ARMSE
ETKFInfl_RMSE = RMSE;
clear RMSE
ETKFInfl_ARelRMSE = ARelRMSE;
clear ARelRMSE
ETKFInfl_RelRMSE = RelRMSE;
clear RelRMSE

save ETKFOptInflTest2

% Opt Inflation: ~1.02-1.03