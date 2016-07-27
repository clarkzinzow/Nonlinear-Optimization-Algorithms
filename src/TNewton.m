function [inform,x] = TNewton(fun,x,nparams)

global numf numg numH;
numf = 0; numg = 0; numH = 0;
n = length(x.p);

% lsparams = struct('ftol',1e-4,'gtol',0.9,'maxit',100);
lsparams = struct('ftol',1e-4,'gtol',0.9,'xtol',1e-6,'stpmin',0,'stpmax',1e20,'maxfev',10000);
[x.f,x.g] = feval(fun,x.p,3);

cgits = 0;
ng = norm(x.g);
if (ng/min(1e3,1 + abs(x.f)) < nparams.toler) 
  inform = struct('iter',0,'status',1,'cgits',cgits);
  return;
end

alfa0 = 1;
s = zeros(n,1);

for iter=1:nparams.maxit
  H = feval(fun,x.p,4);
  [d,its] = cgStep(x.g,H,min(0.5,sqrt(ng))*ng);
  cgits = cgits + its;
  
  [alfa,x] = StepSize(fun, x, d, alfa0, lsparams);
  ng = norm(x.g);
%  fprintf('iter %d f = %g n(g) = %g\n',iter,x.f,ng);
  if (ng/min(1e3,1 + abs(x.f)) < nparams.toler) 
    inform = struct('iter',iter,'status',1,'cgits',cgits);
    return;
  end
end  
inform = struct('iter',iter,'status',0,'cgits',cgits);
return;

function [p,j] = cgStep(g, B, toler)

p = zeros(size(g));
r = g;
d = -r;
nr = norm(r);
if nr < toler
  j = 0;
  return;
end
for j = 1:length(p)
  kap = d'*B*d;
  if (kap <= 0)
    if (j == 1) 
      p = -r;
    end
    return;
  end
  alpha = nr^2/kap;
  p = p + alpha*d;
  r = r + alpha*B*d;
  nrplus = norm(r);
  if nrplus < toler
    return;
  end
  beta = (nrplus/nr)^2;
  nr = nrplus;
  d = -r + beta*d;
end
disp('too many cg iters');
return;
