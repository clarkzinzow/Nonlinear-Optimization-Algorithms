      function [stp,x,inform] = cvsrch(fcn,x,s,stp,params)
%   Translation of minpack subroutine cvsrch
%   Dianne O'Leary   July 1991
%     **********
%
%     Subroutine cvsrch
%
%     The purpose of cvsrch is to find a step which satisfies 
%     a sufficient decrease condition and a curvature condition.
%     The user must provide a subroutine which calculates the
%     function and the gradient.
%
%     At each stage the subroutine updates an interval of
%     uncertainty with endpoints stx and sty. The interval of
%     uncertainty is initially chosen so that it contains a 
%     minimizer of the modified function
%
%          f(x+stp*s) - f(x) - ftol*stp*(gradf(x)'s).
%
%     If a step is obtained for which the modified function 
%     has a nonpositive function value and nonnegative derivative, 
%     then the interval of uncertainty is chosen so that it 
%     contains a minimizer of f(x+stp*s).
%
%     The algorithm is designed to find a step which satisfies 
%     the sufficient decrease condition 
%
%           f(x+stp*s) <= f(x) + ftol*stp*(gradf(x)'s),
%
%     and the curvature condition
%
%           abs(gradf(x+stp*s)'s)) <= gtol*abs(gradf(x)'s).
%
%     If ftol is less than gtol and if, for example, the function
%     is bounded below, then there is always a step which satisfies
%     both conditions. If no step can be found which satisfies both
%     conditions, then the algorithm usually stops when rounding
%     errors prevent further progress. In this case stp only 
%     satisfies the sufficient decrease condition.
%
%     The subroutine statement is
%
%        subroutine cvsrch(fcn,n,x,f,g,s,stp,ftol,gtol,xtol,
%                          stpmin,stpmax,maxfev,info,nfev,wa)
%     where
%
%	fcn is the name of the user-supplied subroutine which
%         calculates the function and the gradient.  fcn must 
%      	  be declared in an external statement in the user 
%         calling program, and should be written as follows.
%
%         function [f,g] = fcn(n,x) (Matlab)     (10/2010 change in documentation)
%	  (derived from Fortran subroutine fcn(n,x,f,g) )
%         integer n
%         f
%         x(n),g(n)
%	  ----------
%         Calculate the function at x and
%         return this value in the variable f.
%         Calculate the gradient at x and
%         return this vector in g.
%	  ----------
%	  return
%	  end
%
%       n is a positive integer input variable set to the number
%	  of variables.
%
%	x is an array of length n. On input it must contain the
%	  base point for the line search. On output it contains 
%         x + stp*s.
%
%	f is a variable. On input it must contain the value of f
%         at x. On output it contains the value of f at x + stp*s.
%
%	g is an array of length n. On input it must contain the
%         gradient of f at x. On output it contains the gradient
%         of f at x + stp*s.
%
%	s is an input array of length n which specifies the
%         search direction.
%
%	stp is a nonnegative variable. On input stp contains an
%         initial estimate of a satisfactory step. On output
%         stp contains the final estimate.
%
%       ftol and gtol are nonnegative input variables. Termination
%         occurs when the sufficient decrease condition and the
%         directional derivative condition are satisfied.
%
%	xtol is a nonnegative input variable. Termination occurs
%         when the relative width of the interval of uncertainty 
%	  is at most xtol.
%
%	stpmin and stpmax are nonnegative input variables which 
%	  specify lower and upper bounds for the step.
%
%	maxfev is a positive integer input variable. Termination
%         occurs when the number of calls to fcn is at least
%         maxfev by the end of an iteration.
%
%	info is an integer output variable set as follows:
%	  
%	  info = 0  Improper input parameters.
%
%	  info = 1  The sufficient decrease condition and the
%                   directional derivative condition hold.
%
%	  info = 2  Relative width of the interval of uncertainty
%		    is at most xtol.
%
%	  info = 3  Number of calls to fcn has reached maxfev.
%
%	  info = 4  The step is at the lower bound stpmin.
%
%	  info = 5  The step is at the upper bound stpmax.
%
%	  info = 6  Rounding errors prevent further progress.
%                   There may not be a step which satisfies the
%                   sufficient decrease and curvature conditions.
%                   Tolerances may be too small.
%
%       nfev is an integer output variable set to the number of
%         calls to fcn.
%
%	wa is a work array of length n.
%
%     Subprograms called
%
%	user-supplied......fcn
%
%	MINPACK-supplied...cstep
%
%	FORTRAN-supplied...abs,max,min
%	  
%     Argonne National Laboratory. MINPACK Project. June 1983
%     Jorge J. More', David J. Thuente
%
%     **********
      p5 = .5;
      p66 = .66;
      xtrapf = 4;
      inform.info = 0;
      infoc = 1;

%
%     Check the input parameters for errors.
%
      if (length(x.p) <= 0 | stp <= 0.0 | params.ftol < 0.0 |  ...
          params.gtol < 0.0 | params.xtol < 0.0 | params.stpmin < 0.0  ...
          | params.stpmax < params.stpmin | params.maxfev <= 0) 
         return
      end
%
%     Compute the initial gradient in the search direction
%     and check that s is a descent direction.
%
      dginit = x.g'*s;
      if (dginit >= 0.0) 
          return
      end
%
%     Initialize local variables.
%
      brackt = 0;
      stage1 = 1;
      inform.nfev = 0;
      finit = x.f;
      dgtest = params.ftol*dginit;
      width = params.stpmax - params.stpmin;
      width1 = 2*width;
      wa = x.p;
%
%     The variables stx, fx, dgx contain the values of the step, 
%     function, and directional derivative at the best step.
%     The variables sty, fy, dgy contain the value of the step,
%     function, and derivative at the other endpoint of
%     the interval of uncertainty.
%     The variables stp, f, dg contain the values of the step,
%     function, and derivative at the current step.
%
      stx = 0.0;
      fx = finit;
      dgx = dginit;
      sty = 0.0;
      fy = finit;
      dgy = dginit;
%
%     Start of iteration.
%
   while (1)   
%
%        Set the minimum and maximum steps to correspond
%        to the present interval of uncertainty.
%
         if (brackt) 
            stmin = min(stx,sty);
            stmax = max(stx,sty);
         else
            stmin = stx;
            stmax = stp + xtrapf*(stp - stx);
         end 
%
%        Force the step to be within the bounds stpmax and stpmin.
%
         stp = max(stp,params.stpmin);
         stp = min(stp,params.stpmax);
%
%        If an unusual termination is to occur then let 
%        stp be the lowest point obtained so far.
%
         if ((brackt & (stp <= stmin | stp >= stmax)) ...
            | inform.nfev >= params.maxfev-1 | infoc == 0 ...
            | (brackt & stmax-stmin <= params.xtol*stmax)) 
            stp = stx;
         end
%
%        Evaluate the function and gradient at stp
%        and compute the directional derivative.
%
         x.p = wa + stp * s;
         [x.f,x.g] = feval(fcn,x.p,3);
         inform.nfev = inform.nfev + 1;
         dg = x.g' * s;
         ftest1 = finit + stp*dgtest;
%
%        Test for convergence.
%
         if ((brackt & (stp <= stmin | stp >= stmax)) | infoc == 0) 
                  inform.info = 6;
         end
         if (stp == params.stpmax & x.f <= ftest1 & dg <= dgtest) 
                  inform.info = 5;
         end
         if (stp == params.stpmin & (x.f > ftest1 | dg >= dgtest)) 
                  inform.info = 4;
         end
         if (inform.nfev >= params.maxfev) 
                  inform.info = 3;
         end
         if (brackt & stmax-stmin <= params.xtol*stmax) 
                  inform.info = 2;
         end
         if (x.f <= ftest1 & abs(dg) <= params.gtol*(-dginit)) 
                  inform.info = 1;
         end
%
%        Check for termination.
%
         if (inform.info ~= 0) 
                  return
         end
%
%        In the first stage we seek a step for which the modified
%        function has a nonpositive value and nonnegative derivative.
%
         if (stage1 & x.f <= ftest1 & dg >= min(params.ftol,params.gtol)*dginit) 
                stage1 = 0;
         end
%
%        A modified function is used to predict the step only if
%        we have not obtained a step for which the modified
%        function has a nonpositive function value and nonnegative 
%        derivative, and if a lower function value has been  
%        obtained but the decrease is not sufficient.
%
         if (stage1 & x.f <= fx & x.f > ftest1) 
%
%           Define the modified function and derivative values.
%
            fm = x.f - stp*dgtest;
            fxm = fx - stx*dgtest;
            fym = fy - sty*dgtest;
            dgm = dg - dgtest;
            dgxm = dgx - dgtest;
            dgym = dgy - dgtest;
% 
%           Call cstep to update the interval of uncertainty 
%           and to compute the new step.
%
            [stx,fxm,dgxm,sty,fym,dgym,stp,fm,dgm,brackt,infoc] ...
             = cstep(stx,fxm,dgxm,sty,fym,dgym,stp,fm,dgm, ...
                     brackt,stmin,stmax);
%
%           Reset the function and gradient values for f.
%
            fx = fxm + stx*dgtest;
            fy = fym + sty*dgtest;
            dgx = dgxm + dgtest;
            dgy = dgym + dgtest;
         else
% 
%           Call cstep to update the interval of uncertainty 
%           and to compute the new step.
%
            [stx,fx,dgx,sty,fy,dgy,stp,x.f,dg,brackt,infoc] ...
             = cstep(stx,fx,dgx,sty,fy,dgy,stp,x.f,dg, ...
                     brackt,stmin,stmax);
         end
%
%        Force a sufficient decrease in the size of the
%        interval of uncertainty.
%
         if (brackt) 
            if (abs(sty-stx) >= p66*width1) 
              stp = stx + p5*(sty - stx);
            end
            width1 = width;
            width = abs(sty-stx);
         end
%
%        End of iteration.
%
     end
%
%     Last card of subroutine cvsrch.
%

     function  [stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt,info] ...
       = cstep(stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt,stpmin,stpmax)
%   Translation of minpack subroutine cstep 
%   Dianne O'Leary   July 1991
%     **********
%
%     Subroutine cstep
%
%     The purpose of cstep is to compute a safeguarded step for
%     a linesearch and to update an interval of uncertainty for
%     a minimizer of the function.
%
%     The parameter stx contains the step with the least function
%     value. The parameter stp contains the current step. It is
%     assumed that the derivative at stx is negative in the
%     direction of the step. If brackt is set true then a
%     minimizer has been bracketed in an interval of uncertainty
%     with endpoints stx and sty.
%
%     The subroutine statement is
%
%       subroutine cstep(stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt,
%                        stpmin,stpmax,info)
% 
%     where
%
%       stx, fx, and dx are variables which specify the step,
%         the function, and the derivative at the best step obtained
%         so far. The derivative must be negative in the direction
%         of the step, that is, dx and stp-stx must have opposite 
%         signs. On output these parameters are updated appropriately.
%
%       sty, fy, and dy are variables which specify the step,
%         the function, and the derivative at the other endpoint of
%         the interval of uncertainty. On output these parameters are 
%         updated appropriately.
%
%       stp, fp, and dp are variables which specify the step,
%         the function, and the derivative at the current step.
%         If brackt is set true then on input stp must be
%         between stx and sty. On output stp is set to the new step.
%
%       brackt is a logical variable which specifies if a minimizer
%         has been bracketed. If the minimizer has not been bracketed
%         then on input brackt must be set false. If the minimizer
%         is bracketed then on output brackt is set true.
%
%       stpmin and stpmax are input variables which specify lower 
%         and upper bounds for the step.
%
%       info is an integer output variable set as follows:
%         If info = 1,2,3,4,5, then the step has been computed
%         according to one of the five cases below. Otherwise
%         info = 0, and this indicates improper input parameters.
%
%     Subprograms called
%
%       FORTRAN-supplied ... abs,max,min,sqrt
%                        ... dble
%
%     Argonne National Laboratory. MINPACK Project. June 1983
%     Jorge J. More', David J. Thuente
%
%     **********
      p66 = 0.66;
      info = 0;
%
%     Check the input parameters for errors.
%
      if ((brackt & (stp <= min(stx,sty) | ...
          stp >= max(stx,sty))) | ...
          dx*(stp-stx) >= 0.0 | stpmax < stpmin) 
         return
      end
%
%     Determine if the derivatives have opposite sign.
%
      sgnd = dp*(dx/abs(dx));
%
%     First case. A higher function value.
%     The minimum is bracketed. If the cubic step is closer
%     to stx than the quadratic step, the cubic step is taken,
%     else the average of the cubic and quadratic steps is taken.
%
      if (fp > fx) 
         info = 1;
         bound = 1;
         theta = 3*(fx - fp)/(stp - stx) + dx + dp;
         s = norm([theta,dx,dp],inf);
         gamma = s*sqrt((theta/s)^2 - (dx/s)*(dp/s));
         if (stp < stx) 
             gamma = -gamma;
         end
         p = (gamma - dx) + theta;
         q = ((gamma - dx) + gamma) + dp;
         r = p/q;
         stpc = stx + r*(stp - stx);
         stpq = stx + ((dx/((fx-fp)/(stp-stx)+dx))/2)*(stp - stx);
         if (abs(stpc-stx) < abs(stpq-stx)) 
            stpf = stpc;
         else
           stpf = stpc + (stpq - stpc)/2;
         end 
         brackt = 1;
%
%     Second case. A lower function value and derivatives of
%     opposite sign. The minimum is bracketed. If the cubic
%     step is closer to stx than the quadratic (secant) step, 
%     the cubic step is taken, else the quadratic step is taken.
%
      elseif (sgnd < 0.0) 
         info = 2;
         bound = 0;
         theta = 3*(fx - fp)/(stp - stx) + dx + dp;
         s = norm([theta,dx,dp],inf);
         gamma = s*sqrt((theta/s)^2 - (dx/s)*(dp/s));
         if (stp > stx) 
            gamma = -gamma;
         end
         p = (gamma - dp) + theta;
         q = ((gamma - dp) + gamma) + dx;
         r = p/q;
         stpc = stp + r*(stx - stp);
         stpq = stp + (dp/(dp-dx))*(stx - stp);
         if (abs(stpc-stp) > abs(stpq-stp))
            stpf = stpc;
         else
            stpf = stpq;
         end 
         brackt = 1;
%
%     Third case. A lower function value, derivatives of the
%     same sign, and the magnitude of the derivative decreases.
%     The cubic step is only used if the cubic tends to infinity 
%     in the direction of the step or if the minimum of the cubic
%     is beyond stp. Otherwise the cubic step is defined to be 
%     either stpmin or stpmax. The quadratic (secant) step is also 
%     computed and if the minimum is bracketed then the the step 
%     closest to stx is taken, else the step farthest away is taken.
%
      elseif (abs(dp) < abs(dx)) 
         info = 3;
         bound = 1;
         theta = 3*(fx - fp)/(stp - stx) + dx + dp;
         s = norm([theta,dx,dp],inf);
%
%        The case gamma = 0 only arises if the cubic does not tend
%        to infinity in the direction of the step.
%
         gamma = s*sqrt(max(0.,(theta/s)^2 - (dx/s)*(dp/s)));
         if (stp > stx) 
             gamma = -gamma;
         end
         p = (gamma - dp) + theta;
         q = (gamma + (dx - dp)) + gamma;
         r = p/q;
         if (r < 0.0 & gamma ~= 0.0)
            stpc = stp + r*(stx - stp);
         elseif (stp > stx)
            stpc = stpmax;
         else
            stpc = stpmin;
         end 
         stpq = stp + (dp/(dp-dx))*(stx - stp);
         if (brackt) 
            if (abs(stp-stpc) < abs(stp-stpq)) 
               stpf = stpc;
            else
               stpf = stpq;
            end
         else
            if (abs(stp-stpc) > abs(stp-stpq)) 
               stpf = stpc;
            else
               stpf = stpq;
            end 
         end 
%
%     Fourth case. A lower function value, derivatives of the
%     same sign, and the magnitude of the derivative does
%     not decrease. If the minimum is not bracketed, the step
%     is either stpmin or stpmax, else the cubic step is taken.
%
      else
         info = 4;
         bound = 0;
         if (brackt) 
            theta = 3*(fp - fy)/(sty - stp) + dy + dp;
            s = norm([theta,dy,dp],inf);
            gamma = s*sqrt((theta/s)^2 - (dy/s)*(dp/s));
            if (stp > sty) 
                gamma = -gamma;
            end
            p = (gamma - dp) + theta;
            q = ((gamma - dp) + gamma) + dy;
            r = p/q;
            stpc = stp + r*(sty - stp);
            stpf = stpc;
         elseif (stp > stx)
            stpf = stpmax;
         else
            stpf = stpmin;
         end 
      end 
%
%     Update the interval of uncertainty. This update does not
%     depend on the new step or the case analysis above.
%
      if (fp > fx) 
         sty = stp;
         fy = fp;
         dy = dp;
      else
         if (sgnd < 0.0)
            sty = stx;
            fy = fx;
            dy = dx;
         end 
         stx = stp;
         fx = fp;
         dx = dp;
      end
%
%     Compute the new step and safeguard it.
%
      stpf = min(stpmax,stpf);
      stpf = max(stpmin,stpf);
      stp = stpf;
      if (brackt & bound)
         if (sty > stx) 
            stp = min(stx+p66*(sty-stx),stp);
         else
            stp = max(stx+p66*(sty-stx),stp);
         end
      end
      return
%
%     Last card of subroutine cstep.
%



