function [inform,x] = SteepDescent(fun,x,sdparams)
%  SteepDescent implements steepest descent using
%  the subroutine StepSize.m.
%
%  Implementation Parameters:
%    * STOP when
%      - Objective function has norm less than gtol.
%      - maxit steps have been made.
%    * ftol = 1.0e-20
%    * gtol = 1.0e-4
%    * xtol = 1.0e-20
%    * maxit = 1000
%    * p = -grad(f(x_k))
%    * alfa = max(10*params.xtol,
%             f'(x_{k-1}:d_{k-1})*alfa/f'(x_k:d_k))
%    * geoparams = parameter structure supplied if using geodesic
%                  function
%
%  Input:
%    fun      - a pointer to a function
%    x        - the following structure:
%               * x.p - the starting point values
%               * x.f - the function value of x.p
%               * x.g - the gradient value of x.p
%    sdparams - the following structure, as an example:
%         sdparams = struct('maxit',1000,'toler',1.0e-4,
%                           'geoparams',geoparams);
%
%  Output:
%    inform - structure containing two fields:
%      * inform.status - 1 if gradient tolerance was
%                        achieved;
%                        0 if not.
%      * inform.iter   - the number of steps taken
%    x      - the solution structure, with point, function,
%             and gradient evaluations at the solution point.

global numf numg

numf = 0;
numg = 0;

params = struct('ftol', 1.0e-20, 'gtol', sdparams.toler,...
                'xtol', 1.0e-20);  % Parameter structure

% If calling geodesic function...
if isfield(sdparams, 'geoparams') && isstruct(sdparams.geoparams)
    params.geoparams = sdparams.geoparams;  % Add geoparams to params
    % alfa = sdparams.geoparams.alpha;
    x.f = feval(fun, x.p, 1, params.geoparams);  % Function value at starting point
    x.g = feval(fun, x.p, 2, params.geoparams);  % Gradient value at starting point
else
    params.geoparams = -1;  % If not, flag geoparams
    x.f = feval(fun, x.p, 1);  % Function value at starting point
    x.g = feval(fun, x.p, 2);  % Gradient value at starting point
end

alfa = 1;  % Initial alpha value
inform = struct('status', 0, 'iter', 0);  % Information structure
while inform.iter < sdparams.maxit && norm(x.g) >= sdparams.toler
    pg = x.g;  % Store previous gradient value
    d = -1*x.g;  % Steepest descent direction
    [alfa,x] = StepSize(fun, x, d, alfa, params);  % Get step size
    alfa = max(10*params.xtol, pg'*d*alfa/(x.g'*(-1*x.g)));  % New alpha
    inform.iter = inform.iter + 1;  % Increment step counter
end
if norm(x.g) < sdparams.toler  % If gradient norm less than tolerance
    inform.status = 1;  % Then success
end
return;  % Return inform and final point x
end