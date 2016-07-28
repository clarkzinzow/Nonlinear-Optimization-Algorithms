function [inform, x] = SteepDescent(fun, x, sdparams)
%  Implements steepest descent using simple Wolfe conditions via StepSize.m.
%
%  Implementation Parameters:
%    * STOP when
%      - Objective function gradient has norm less than gtol.
%      - maxit steps have been taken.
%    * ftol = 1.0e-20
%    * gtol = 1.0e-4
%    * xtol = 1.0e-20
%    * maxit = 1000
%    * p_k = -grad(f(x_k))
%    * alfa_k = max(10*xtol, f'(x_{k-1})*p_{k-1}*alfa_{k-1}/(f'(x_k)*p_k))
%
%  Input:
%    fun      - a pointer to a function
%    x        - the following structure:
%               * x.p - the starting point values
%               * x.f - the function value of x.p
%               * x.g - the gradient value of x.p
%    sdparams - the following structure, as an example:
%         sdparams = struct('maxit',1000,'toler',1.0e-4);
%
%  Output:
%    inform - structure containing two fields:
%      * inform.status - 1 if gradient tolerance was
%                        achieved;
%                        0 if not.
%      * inform.iter   - the number of steps taken
%    x      - the solution structure, with point, function,
%             and gradient evaluations at the solution point.

%  Number of function and gradient evaluations.
global numf numg
numf = 0;
numg = 0;

%  Populate local caching of nparams parameters, and params parameters.
toler = sdparams.toler;  % Set gradient tolerance.
maxit = sdparams.maxit;  % Set maximum number of allowed iterations.
xtol = 1.0e-20;  % Set point tolerance.
ftol = 1.0e-20;  % Set function tolerance.
gtol = toler;  % Set gradient tolerance.

%  Initialize parameter structure for StepSize function call.
params = struct('ftol', ftol, 'gtol', gtol, 'xtol', xtol);

% Below was tailored to geodesic functions.
% If calling geodesic function...
% if isfield(sdparams, 'geoparams') && isstruct(sdparams.geoparams)
%     params.geoparams = sdparams.geoparams;
%     x.f = feval(fun, x.p, 1, params.geoparams);
%     x.g = feval(fun, x.p, 2, params.geoparams);
% end
% params.geoparams = -1;

%  Compute function and gradient at starting point.
x.f = feval(fun, x.p, 1);
x.g = feval(fun, x.p, 2);

alfa = 1;  % Initial alpha value.
iter = 0;  % Number of iterations.
while iter < maxit && norm(x.g) >= toler
    pg = x.g;  % Store previous gradient value.
    d = -x.g;  % Steepest descent direction.
    
    %  Get step size satisfying simple Wolfe conditions.
    [alfa, x] = StepSize(fun, x, d, alfa, params);
    alfa = max(10*xtol, pg'*d*alfa/(x.g'*(-x.g)));  % Calculate new alpha.
    iter = iter + 1;  % Increment step counter.
end
inform.iter = iter;  % Number of iterations = iter at this point.

%  If norm of gradient is less than tolerance, the method succeeded.
if norm(x.g) < toler
    inform.status = 1;  % Update status to success indicator, 1.
else
    inform.status = 0;  % Update status to failure indicator, 0.
end
return;  % Return inform and final point x
end