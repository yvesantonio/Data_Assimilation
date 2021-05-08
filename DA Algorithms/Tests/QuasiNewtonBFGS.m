function [xk, fk, counter, error, xks, fks] = ...
    QuasiNewtonBFGS(f, x0, epsilon, maxiterations)

xks = [x0'];
[fks, dFk, B] = feval(f,x0);

xk = x0;
counter = 0;
error = 1e300;

Hk = B^-1;
I = eye(length(x0));
while error > epsilon && counter < maxiterations
    pk = -Hk*dFk;
    alpha = fminsearch(@(a) feval(f, xk + a*pk), 0.0); %exact 1-d opt
    
    sk = alpha*pk;
    xk = xk + sk;
    [fx, dFk1] = feval(f, xk);
    yk = dFk1 - dFk;
    rhok = 1/(yk'*sk);
    Hk = (I - rhok*sk*yk')*Hk*(I - rhok*yk*sk') + rhok*(sk*sk');
    dFk = dFk1;
    
    error = norm(dFk);
    fk = feval(f,xk);
    
    xks = [xks; xk'];
    fks = [fks; fk];
    counter = counter + 1;
end

end %of function
