function [xk, fk, counter, error, xks, fks] = ...
    FletcherReeves(f, x0, epsilon, maxiterations)

xks = [x0'];
[fks, dFks] = feval(f,x0);

xk = x0;
counter = 0;
error = 1e300;

dF1 = -dFks;
pk  = -dFks;
while error > epsilon && counter < maxiterations
    
    counter = counter + 1;
    alpha = fminsearch(@(a) feval(f, xk + a*pk), 1); %exact 1-d opt
    alpha 
    xk = xk + alpha*pk;
    dF = dF1;
    [fx, dF1] = feval(f, xk);
    betafr = (dF1'*dF1)/(dF'*dF);
    pk = -dF1 + betafr*pk;
    
    error = norm(dF1);
    fk = feval(f,xk);
    
    xks = [xks; xk'];
    fks = [fks; fk];
    
end

end %of function
