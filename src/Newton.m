function [inform, x] = Newton(func, x, nparams)
%  Implements Newton's method, which uses the Newton iteration for the search
%  direction:
%
%                      p = -(Hessian(f)^(-1)) grad(f).
%
%  Moreover, we either use a 'direct' modified Hessian method (if 'direct' is
%  supplied as the 'method' in nparams) or we use a Cholesky with added
%  multiple of the identity method for Hessian modification.  We are therefore
%  able to deal with functions whose Hessian may not be positive definite away
%  from the solution (entailing a non-descent direction for the Newton
%  direction.)
%
%  Input:
%    func      - a pointer to a function
%    x        - the following structure:
%               * x.p - the starting point values
%               * x.f - the function value of x.p
%               * x.g - the gradient value of x.p
%               * x.h - the Hessian value of x.p
%               with only x.p guaranteed to be set.
%    nparams - the following structure, as an example:
%         nparams = struct('maxit',1000,'toler',1.0e-4,'method','direct');
%
%  Output:
%    inform - structure containing two fields:
%      * inform.status - 1 if gradient tolerance was
%                        achieved;
%                        0 if not.
%      * inform.iter   - the number of steps taken
%    x      - the solution structure, with the solution point along with
%             function and gradient evaluations thereof.

%  Number of function, gradient, and Hessian evaluations.
global numf numg numh
numf = 0;
numg = 0;
numh = 0;

%  Populate local caching of nparams parameters.
toler = nparams.toler;  % Set gradient tolerance.
maxit = nparams.maxit;  % Set maximum number of allowed iterations.
method = nparams.method;  % Set method.
xc.p = x.p;  % Set the current point to the initial point, x.p.

%  Initialize parameter structure for StepSize function call.
params = struct('ftol', 1e-4, 'gtol', 0.9, 'xtol', 1e-6, 'stpmin', 0, ...
                'stpmax', 1e20, 'maxfev', 10000);

for i = 1:maxit
    %  Compute function, gradient, and Hessian at current point.
    xc.f = feval(func, xc.p, 1);
    xc.g = feval(func, xc.p, 2);
    xc.h = sparse(feval(func, xc.p, 4));
    
    % Check for the termination condition: norm of gradient less than toler.
    if norm(xc.g) < toler
        inform.status = 1;  % Indicates success.
        inform.iter = i;  % Number of iterations.
        x.p = xc.p;
        x.f = xc.f;
        x.g = xc.g;
        return;
    end
    
    % Use the direct method.
    if strcmp(method, 'direct')
        s = -xc.h \ xc.g;  % Search direction.
        
        % If this is not a descent direction...
        if xc.g' * s >= 0
            D = zeros(size(xc.h));  % Modified Hessian matrix.
            for j = 1 : size(xc.p, 1)
                D(j,j) = 1.0 / max(0.01, abs(xc.h(j,j)));
            end
            %  D is now a positive diagonal matrix with a diagonal scaled
            %  by the inverse of the corresponding diagonal element of the
            %  Hessian (or by 100 if 1/abs(xc.h(j,j)) > 100.)
            s = -D*xc.g;  % New search direction that is a descent direction.
        end
    else
        % Hessian modification method: Cholesky with added multiple of identity.
        [R, ~] = CholeskyMultIdentity(xc.h);
        s = -R \ (R'\xc.g);  % Search direction.
    end
    
    %  Get step size that satisfies simple Wolfe conditions.
    [alfa, x] = StepSize(func, xc, s, 1, params);
    %  Update current point in s-direction.
    xc.p = xc.p + alfa * s;
end
%  If reached, method failed.
inform.status = 0;  % Update status to failure indicator, 0.
inform.iter = maxit;  % Number of iterations = i = maxit at this point.
x.p = xc.p;
x.f = feval(func, x.p, 1);
x.g = feval(func, x.p, 2);
return;  % Return inform and final point.
end


