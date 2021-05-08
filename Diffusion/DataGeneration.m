clear
close all

load NominalU.mat

%B = diag([0.1*ones(1, 15), ones(1, 20) , 0.1*ones(1, 15)]);
B = 5*eye(nx);

u_0_bg = (mvnrnd(u_0_t, B))';
u_0_bg(1) = 1;
u_0_bg(end) = 1;

u_0_en = ensembleInit(u_0_bg, B, 20);

u_meas = zeros(nx-2, nt);
R = 0.1*eye(nx-2);

for i=2:nt+1
    u_meas(:,i-1) = (mvnrnd(u_tot(2:nx-1,i), R))';
end

M = @(t, x) [UL; D\(x(2:nx-1)+bc); UR] ;

u_bg_tot = zeros(nx, nt+1);
u_bg_tot(:,1) = u_0_bg;

for i=2:nt+1
    u_bg_tot(:,i) = M(i, u_bg_tot(:,i-1));    
end

ArelRMSE_bg = averageRelativeRootMeanSquareError(u_bg_tot, u_tot)
ArelRMSE_meas = averageRelativeRootMeanSquareError(u_meas, u_tot(2:nx-1,2:end))

H = @(t, x) x(2:nx-1);

selector = cell(1,nx);

parfor i1=1:nx
    selector{i1} = simpleSelector(i1, nx-2, 3);
end
sel = @(i) selector{i};

save PerturbedU