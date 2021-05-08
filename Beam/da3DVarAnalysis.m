clear;
close all;

load 3DVarData

state_meas = state_meas';

ARelRMSE_bg = averageRelativeRootMeanSquareError(state_bg, state)
ARelRMSE_meas = averageRelativeRootMeanSquareError(state_meas, state)

n = length(state_bg);
H = eye(n);

% % % [ state_a, J_xa] = da_var_3DVar(state_bg, state_meas, eye(n), R, H, @(f, x0) FletcherReeves(f, x0, 1e-3, 250));
% % % ARelRMSE_bg = averageRelativeRootMeanSquareError(state_bg, state)
% % % ARelRMSE_meas = averageRelativeRootMeanSquareError(state_meas, state)
% % % ARelRMSE_a = averageRelativeRootMeanSquareError(state_a, state)

% % % % % % % % ARMSE_bg = averageRootMeanSquareError(state_bg, state)
% % % % % % % % ARMSE_meas = averageRootMeanSquareError(state_meas, state)
% % % % % % % % ARMSE_a = averageRootMeanSquareError(state_a, state)

state_a = state_bg + B*H'*(H*B*H' + R)^-1*(state_meas - H*state_bg);
ARelRMSE_bg = averageRelativeRootMeanSquareError(state_bg, state)
ARelRMSE_meas = averageRelativeRootMeanSquareError(state_meas, state)
ARelRMSE_a = averageRelativeRootMeanSquareError(state_a, state)

save 3DVarAnalitic