function [ rho ] = simpleRho(n, c)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

rho = diag(ones(n,1));    
for i=1:c-1
    r = i/c;
    Gr = 1 - 5/3*r^2 + 5/8*r^3 + r^4/2 - r^5/4;
    
    for j = 1:i        
        rho(j, n-i+j) = Gr;
        rho(n-i+j, j) = Gr;        
    end
    
    Gr = Gr.*ones(n-i,1);
    rho = rho + diag(Gr, i);
    rho = rho + diag(Gr, -i);
end

for i=c:((c*2)-1)
    r = i/c;
    Gr = 4 - 5*r + 5/3*r^2 + 5/8*r^3- r^4/2 + r^5/12 - 2/(3*r);
    
	for j = 1:i        
        rho(j, n-i+j) = Gr;
        rho(n-i+j, j) = Gr;        
    end
    
    Gr = Gr.*ones(n-i,1);
    rho = rho + diag(Gr, i);
    rho = rho + diag(Gr, -i);    
end
end
