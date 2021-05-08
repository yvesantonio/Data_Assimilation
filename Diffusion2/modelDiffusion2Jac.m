function [ J_out ] = modelDiffusion2Jac(nx, ny, dx, dy, dt)

vis=0.1;                         %Diffusion coefficient/viscocity
%Calculating the coefficient matrix for the implicit scheme
Ex=sparse(2:nx-2,1:nx-3,1,nx-2,nx-2);
Ax=Ex+Ex'-2*speye(nx-2);        %Dirichlet B.Cs
%Ax(1,1)=-1; Ax(nx-2,nx-2)=-1;  %Neumann B.Cs
Ey=sparse(2:ny-2,1:ny-3,1,ny-2,ny-2);
Ay=Ey+Ey'-2*speye(ny-2);        %Dirichlet B.Cs
%Ay(1,1)=-1; Ay(ny-2,ny-2)=-1;  %Neumann B.Cs
A=kron(Ay/dy^2,speye(nx-2))+kron(speye(ny-2),Ax/dx^2);
D=speye((nx-2)*(ny-2))-vis*dt*A;

J_out = sparse([ zeros(nx+ny-2, nx*ny); zeros((nx-2)*(ny-2), nx+ny-2) D^-1 zeros((nx-2)*(ny-2), nx+ny-2); zeros(nx+ny-2, nx*ny) ]);