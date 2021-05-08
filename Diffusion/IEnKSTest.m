clear
close all

load PerturbedU


L = 3;

u_a_tot = da_seq_bundleIterativeEnsembleKalmanSmoother(u_0_en, u_meas, M, H, R, bIEnKSOptions('OptAlg', @(f, x0) PolakRibiere(f, x0, 1e-6, nx), 'Inflation', 1.01, 'L', L, 'S', 1, 'MDA', 1/(L+1)*ones(L+1,1)));
% u_a_tot = da_seq_bundleIterativeEnsembleKalmanSmoother(u_0_en, u_meas, M, H, R, bIEnKSOptions('OptAlg', @(f, x0) PolakRibiere(f, x0, 1e-6, nx), 'Inflation', 1.01, 'L', L, 'S', 1));
% u_a_tot = da_seq_bundleIterativeEnsembleKalmanSmoother(u_0_en, u_meas, M, H, R, bIEnKSOptions('OptAlg', @(f, x0) PolakRibiere(f, x0, 1e-6, nx), 'Inflation', 1.01, 'Localization', sel, 'L', L, 'S', 1));

ArelRMSE_a = averageRelativeRootMeanSquareError(u_a_tot, u_tot)

save IEnKSMDATest