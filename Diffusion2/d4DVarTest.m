clear
close all

load PerturbedU

J_M = @(t, x) modelDiffusion2Jac(nx, ny, dx, dy, dt);

[ u_a, Jx_a] = da_var_4DVarStrong(u_0_bg, u_meas, M, J_M, H, J_H, B, R, @(f, x0) PolakRibiere(f, x0, 1e-6, 10));

u_a_tot = zeros(nx*ny, nt+1);
u_a_tot(:,1) = u_a;

for i=2:nt+1    
    u_a_tot(:,i) = M(i, u_a_tot(:,i-1));
end

ArelRMSE_a = averageRelativeRootMeanSquareError(u_a_tot, u_tot)

save 4DVarTest