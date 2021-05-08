clear
close all

load PerturbedU


L = 2;

% u_a_tot = da_seq_bundleIterativeEnsembleKalmanSmoother(u_0_en, u_meas, M, H, R, bIEnKSOptions('OptAlg', @(f, x0) PolakRibiere(f, x0, 1e-6, 10), 'Inflation', 1.01, 'L', L, 'S', 1, 'MDA', 1/(L+1)*ones(L+1,1)));
u_a_tot = da_seq_bundleIterativeEnsembleKalmanSmoother(u_0_en, u_meas, M, H, R, bIEnKSOptions('OptAlg', @(f, x0) PolakRibiere(f, x0, 1e-6, 10), 'L', L, 'S', 1, 'Localization', sel));


ArelRMSE_a = averageRelativeRootMeanSquareError(u_a_tot, u_tot)

save IEnKSLocTestL2

clear
close all

load PerturbedU


L = 5;

% u_a_tot = da_seq_bundleIterativeEnsembleKalmanSmoother(u_0_en, u_meas, M, H, R, bIEnKSOptions('OptAlg', @(f, x0) PolakRibiere(f, x0, 1e-6, 10), 'Inflation', 1.01, 'L', L, 'S', 1, 'MDA', 1/(L+1)*ones(L+1,1)));
u_a_tot = da_seq_bundleIterativeEnsembleKalmanSmoother(u_0_en, u_meas, M, H, R, bIEnKSOptions('OptAlg', @(f, x0) PolakRibiere(f, x0, 1e-6, 10), 'L', L, 'S', 1, 'Localization', sel));


ArelRMSE_a = averageRelativeRootMeanSquareError(u_a_tot, u_tot)

save IEnKSLocTestL5