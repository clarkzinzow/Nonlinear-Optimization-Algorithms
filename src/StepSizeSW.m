function [alpha,xp] = StepSizeSW(f,x,d,alpha,params)
% function [alpha,xp] = StepSizeSW(f,x,d,alpha,params)
% Line search algorithm satisfying strong Wolfe conditions
% Algorithms 3.5 on pages 60-61 in Nocedal and Wright
% Requires x.p, x.f and x.g to be initialized

alpha0 = params.stpmin;
c1 = params.ftol;
c2 = params.gtol;
% alpha is alpha_i
gxd = x.g'*d;
% alphap is alpha_{i-1}
alphap = alpha0;
fxp = x.f;
i=1;
while norm(alphap-alpha) > params.xtol
  xp.p = x.p + alpha*d;
  xp.f = feval(f,xp.p,1);
  if (xp.f > x.f + c1*alpha*gxd) | ((i > 1) & (xp.f >= fxp)),
    [alpha,xp] = zoom(f,x,d,alphap,alpha,fxp,c1,c2);
    return;
  end
  xp.g = feval(f,xp.p,2); gxpd = xp.g'*d;
  if abs(gxpd) <= -c2*gxd,
    return;
  end
  if gxpd >= 0,
    [alpha,xp] = zoom(f,x,d,alpha,alphap,xp.f,c1,c2);
    return;
  end
  alphap = alpha;
  fxp = xp.f;
  %  alpha = alpha + (params.stpmax-alpha)*rand(1);
  alpha = alpha + (params.stpmax-alpha)*0.5;
  i = i+1;
end
alpha,
error('No stepsize found');

function [alpha,xp] = zoom(f,x,d,alphal,alphah,fxl,c1,c2)
% function [alpha,xp] = zoom(f,x,d,alphal,alphah,fxl)
% Algorithm 3.6 on page 61 in Nocedal and Wright

gxd = x.g'*d;

while 1
   alpha = 1/2*(alphal+alphah);
   xp.p = x.p + alpha*d;
   xp.f = feval(f,xp.p,1); 
   if ((xp.f > x.f + c1*alpha*gxd) | (xp.f >= fxl)),
      alphah = alpha;
   else
      xp.g = feval(f,xp.p,2); gxpd = xp.g'*d;
      if abs(gxpd) <= -c2*gxd,
        return;
      end
      if gxpd*(alphah-alphal) >= 0,
        alphah = alphal;
      end
      alphal = alpha;
      fxl = xp.f;
   end
end 