% Stub of a descent method
function [xk, fk, counter, error, xks, fks] = ...
    descentmethod(f, x0, epsilon, maxiterations)

  xks = [x0'];
  fks = [feval(f,x0)];

  xk = x0;
  counter = 0;
  error = 1e300;
  
  while error > epsilon && counter < maxiterations

    counter = counter + 1;
    d = -DJ(xk);
    alpha = fminsearch(@(a) feval(f,xk + a*d), 0.0); %exact 1-d opt
    xk = xk + alpha*d;

    error = norm(DJ(xk));
    fk = feval(f,xk');
    
    xks = [xks; xk'];
    fks = [fks; fk];
    
  end

end %of function
