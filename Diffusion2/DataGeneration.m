clear
close all

load NominalU.mat

% % % % % M = @(t, x) modelDiffusion2(t, x, nx, ny, dx, dy, dt);
% % % % % u_test = zeros(2000,31);
% % % % % u_test(:,1) = u_0_t;
% % % % % for i=1:nt
% % % % %     u_test(:,i+1) = M(i+1, u_test(:,i));
% % % % % end
% % % % % ArelRMSE_test = averageRelativeRootMeanSquareError(u_test, u_tot)

B = eye(nx*ny);
u_0_t = u_tot(:,1);
u_0_bg = (mvnrnd(u_0_t, B))';
for i= 1
    U_0_bg = reshape(u_0_bg, nx, ny);
    
    H = ones(nx, ny);
    H (1,:) = 0;
    H (nx,:) = 0;
    H (:,1) = 0;
    H (:,ny) = 0;
    H = reshape(H, [], 1);
    H = diag(H);
    
    j = 1;
    while j <= length(H(:,1)) 
        if(sum(H(j,:)) == 0)
            H = [ H(1:j-1,:); H(j+1:end,:)];
        else
            j = j+1;
        end        
    end
    J_H = H;
    H = @(t, x) J_H*x;
    J_H = @(t, x) J_H;
    
    U_0_bg (1,:)=UW;
    U_0_bg (nx,:)=UE;
    U_0_bg (:,1)=US;
    U_0_bg (:,ny)=UN;
    
    u_0_bg = reshape(U_0_bg, [], 1);
end
u_0_en = ensembleInit(u_0_bg, B, 500);

M = @(t, x) modelDiffusion2(t, x, nx, ny, dx, dy, dt);

u_meas = zeros((nx-2)*(ny-2), nt);
R = 0.005*eye((nx-2)*(ny-2));

u_interior_vector = zeros((nx-2)*(ny-2),nt);
for i=2:nt+1
    U = reshape(u_tot(:,i), nx, ny);
    u_interior_vector(:,i) = reshape(U(2:nx-1,2:ny-1), [], 1);
    u_meas(:,i-1) = (mvnrnd(u_interior_vector(:,i), R))';
end

u_bg_tot = zeros(nx*ny, nt+1);
u_bg_tot(:,1) = u_0_bg;

for i=2:nt+1
    u_bg_tot(:,i) = M(i, u_bg_tot(:,i-1));    
end

ArelRMSE_bg = averageRelativeRootMeanSquareError(u_bg_tot, u_tot)
ArelRMSE_meas = averageRelativeRootMeanSquareError(u_meas, u_interior_vector(:,2:end))

selector = cell(1,nx*ny);

parfor i1=1:nx*ny
    selector{i1} = simpleSelector(i1, (nx-2)*(ny-2), 200);
    selector{i1} = sparse(selector{i1});
end
sel = @(i) selector{i};


save PerturbedU