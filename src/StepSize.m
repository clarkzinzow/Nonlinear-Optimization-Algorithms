function [alfa,x] = StepSize(fun, x, d, alfa, params)
%  Implements simple Wolfe conditions.

x0 = x.p;
Dphi0 = x.g'*d;
if ( (alfa <= 0) || (Dphi0 > 0) )
  error('Initialization of step incorrect');
end;

c1 = params.ftol;
c2 = params.gtol;

phi0 = x.f;
alfaL = 0;
alfaR = inf;

iter = 0;
while abs(alfaR-alfaL) > params.xtol
  iter = iter + 1;
  x.p = x0 + alfa*d;
   
  x.f = feval(fun,x.p,1);
  
  if (x.f >= phi0 + alfa*c1*Dphi0) 
    alfaR = alfa;
    alfa = Interp(alfaL, alfaR);
  else
    x.g = feval(fun,x.p,2);
    DphiAlfa = x.g'*d;
    
    if (DphiAlfa >= c2*Dphi0)
        return;
    else
      alfaL = alfa;
    end
    if isinf(alfaR)
      alfa = Extrap(alfaL);
    else
      alfa = Interp( alfaL, alfaR);
    end	
  end
end
fprintf('step size criteria were not met\n');
fprintf('after %d step size iterations.\n', iter);
error('STEP SIZE FAILURE');
return;

function mid = Interp(left, right)
mid = (left + right)/2.0;
return;

function mid = Extrap(left)
mid = 2*left;
return;

