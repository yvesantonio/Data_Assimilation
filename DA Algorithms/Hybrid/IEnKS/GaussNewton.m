function [xk, fk, counter, error, xks, fks] = ...
    GaussNewton(f, x0, epsilon, maxiterations)

xks = [x0'];
[fks, dFk, Hk] = feval(f,x0);

xk = x0;
counter = 0;
error = 1e300;

while error > epsilon && counter < maxiterations
    dk = Hk\dFk;    
        
    xk = xk - dk;
    [fk, dFk, Hk] = feval(f, xk);
    
    error = norm(dk);    
    
    xks = [xks; xk'];
    fks = [fks; fk];    
    counter = counter + 1;
end

end %of function
