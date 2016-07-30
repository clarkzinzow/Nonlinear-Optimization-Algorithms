function [inform, x] = BFGS(func, x, bfgsparams)
%  Implements the Broyden-Fletcher-Goldfarb-Shanno method, an iterative
%  quasi-Newton algorithm for solving unconstrained nonlinear optimization
%  problems.
%
%  Input:
%    func      - a pointer to a function
%    x        - the following structure:
%               * x.p - the starting point values
%               * x.f - the function value of x.p
%               * x.g - the gradient value of x.p
%               * x.h - the Hessian value of x.p
%               with only x.p guaranteed to be set.
%    bfgsparams - the following structure, as an example:
%         bfgsparams = struct('maxit',1000,'toler',1.0e-4,);
%
%  Output:
%    inform - structure containing two fields:
%      * inform.status - 1 if gradient tolerance was
%                        achieved;
%                        0 if not.
%      * inform.iter   - the number of steps taken
%    x      - the solution structure, with the solution point along with the
%             function evaluation thereof.

%  Number of function, gradient, and Hessian evaluations.
global numf numg numh
numf = 0;
numg = 0;
numh = 0;

%  Populate local caching of bfgsparams parameters.
toler = bfgsparams.toler;  % Set gradient tolerance.
maxit = bfgsparams.maxit;  % Set maximum number of allowed iterations.
xc.p = x.p;  % Set the current point to the initial point, x.p.

%  Initialize parameter structure for StepSize function call.
params = struct('ftol', 1e-4, 'gtol', 0.9, 'xtol', 1e-6, 'stpmin', 0, ...
                'stpmax', 1e20, 'maxfev', 10000);

I = eye(size(xc.p, 1));  % Locally stored identity matrix.
for i = 1:maxit
    %  Compute function and gradient at current point.
    xc.f = feval(func, xc.p, 1);
    xc.g = feval(func, xc.p, 2);
    
    %  Check for termination condition: (scaled) norm of gradient less
    %  than toler.
    if norm(xc.g) / min(1000, 1 + abs(xc.f)) < toler
        inform.status = 1;  % Indicates success.
        inform.iter = i;  % Number of iterations.
        x.p = xc.p;
        x.f = xc.f;
        return;
    end
    
    %  For the first step, we use the identity matrix as an initial inverse
    %  Hessian approximation.
    if i == 1
        H = I;
    else
        %  Update the current inverse Hessian approximation, H.
        H = (I - rho*s*y') * H * (I - rho*y*s') + rho*s*s';
    end
    
    %  Compute the current search direction.
    p = -H * xc.g;
    
    %  Get step size that satisfies simple Wolfe conditions.
    %  NOTE:  alfa = 1 should always be tried first since this step length will
    %         eventually always be accepted (under certain conditions), thereby
    %         producing superlinear convergence of the overall algorithm.
    %         See page 142 of Nocedal and Wright.
    [alfa, x] = StepSize(func, xc, p, 1, params);
    %  Update current point in p-direction with step size alpha.
    xc.p = xc.p + alfa * p;

    %  Update parameters.
    s = alfa * p;  %  s = new_point - prev_point = alfa * p
    y = feval(func, xc.p, 2) - xc.g;  % y = grad(new_point) - grad(prev_point)
    rho = 1 / (y'*s);
    
    %  If in first step, apply inverse Hessian approximation heuristic given by
    %  (6.20) on page 143 of Nocedal and Wright.
    if i == 1
        H = (s'*y)/(y'*y) * I;
    end
end
%  If reached, method failed.
inform.status = 0;  % Update status to failure indicator, 0.
inform.iter = maxit;  % Number of iterations.
x.p = xc.p;
x.f = feval(func, x.p, 1);
x.g = feval(func, x.p, 2);
return;  % Return inform and final point x
end
