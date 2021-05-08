function [ u ] = modelDiffusion2(t, u, nx, ny, dx, dy, dt)

u=reshape(u, nx, ny);            
un=zeros(nx,ny);                 %Preallocating un
vis=0.1;                         %Diffusion coefficient/viscocity
UW=0;                            %x=0 Dirichlet B.C 
UE=0;                            %x=L Dirichlet B.C 
US=0;                            %y=0 Dirichlet B.C 
UN=0;                            %y=L Dirichlet B.C 
UnW=0;                           %x=0 Neumann B.C (du/dn=UnW)
UnE=0;                           %x=L Neumann B.C (du/dn=UnE)
UnS=0;                           %y=0 Neumann B.C (du/dn=UnS)
UnN=0;                           %y=L Neumann B.C (du/dn=UnN)

%B.C vector
bc=zeros(nx-2,ny-2);
bc(1,:)=UW/dx^2; bc(nx-2,:)=UE/dx^2;  %Dirichlet B.Cs
bc(:,1)=US/dy^2; bc(:,ny-2)=UN/dy^2;  %Dirichlet B.Cs
%bc(1,:)=-UnW/dx; bc(nx-2,:)=UnE/dx;  %Neumann B.Cs
%bc(:,1)=-UnS/dy; bc(:,nx-2)=UnN/dy;  %Neumann B.Cs
%B.Cs at the corners:
bc(1,1)=UW/dx^2+US/dy^2; bc(nx-2,1)=UE/dx^2+US/dy^2;
bc(1,ny-2)=UW/dx^2+UN/dy^2; bc(nx-2,ny-2)=UE/dx^2+UN/dy^2;
bc=vis*dt*bc;

%Calculating the coefficient matrix for the implicit scheme
Ex=sparse(2:nx-2,1:nx-3,1,nx-2,nx-2);
Ax=Ex+Ex'-2*speye(nx-2);        %Dirichlet B.Cs
%Ax(1,1)=-1; Ax(nx-2,nx-2)=-1;  %Neumann B.Cs
Ey=sparse(2:ny-2,1:ny-3,1,ny-2,ny-2);
Ay=Ey+Ey'-2*speye(ny-2);        %Dirichlet B.Cs
%Ay(1,1)=-1; Ay(ny-2,ny-2)=-1;  %Neumann B.Cs
A=kron(Ay/dy^2,speye(nx-2))+kron(speye(ny-2),Ax/dx^2);
D=speye((nx-2)*(ny-2))-vis*dt*A;

un=u;

%Implicit method:
U=un;U(1,:)=[];U(end,:)=[];U(:,1)=[];U(:,end)=[];
U=reshape(U+bc,[],1);
U=D\U;
U=reshape(U,nx-2,ny-2);
u(2:nx-1,2:ny-1)=U;
%Boundary conditions
%Dirichlet:
u(1,:)=UW;
u(nx,:)=UE;
u(:,1)=US;
u(:,ny)=UN;

u = reshape(u, [], 1);