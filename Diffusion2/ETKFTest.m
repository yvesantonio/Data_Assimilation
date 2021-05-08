clear
close all

load PerturbedU

u_a_tot = da_seq_EnsembleTransformKalmanFilter(u_0_en, u_meas, M, H, R, 1:nt+1, 'Localization', sel);

ArelRMSE_a = averageRelativeRootMeanSquareError(u_a_tot, u_tot)

save ETKFLocTest