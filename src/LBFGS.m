function [inform, x] = LBFGS(func, x, lbfgsparams)
%  Implements the limited-memory Broyden-Fletcher-Goldfarb-Shanno algorithm, a
%  quasi-Newton method that approximates the Broyden-Fletcher-Goldfarb-Shanno
%  method using a limited amount of memory.

%  Input:
%    func      - a pointer to a function
%    x        - the following structure:
%               * x.p - the starting point values
%               * x.f - the function value of x.p
%               * x.g - the gradient value of x.p
%               * x.h - the Hessian value of x.p
%               with only x.p guaranteed to be set.
%    lbfgsparams - the following structure, as an example:
%         lbfgsparams = struct('maxit', 1000, 'toler', 1.0e-4, 'm', 3);
%
%  Output:
%    inform - structure containing two fields:
%      * inform.status - 1 if gradient tolerance was
%                        achieved;
%                        0 if not.
%      * inform.iter   - the number of steps taken
%    x      - the solution structure, with the solution point along with
%             function, gradient, and Hessian evaluations thereof.

%  Number of function, gradient, and Hessian evaluations.
global numf numg numh
numf = 0;
numg = 0;
numh = 0;

%  Populate local caching of bfgsparams parameters.
toler = lbfgsparams.toler;  % Set gradient tolerance.
maxit = lbfgsparams.maxit;  % Set maximum number of allowed iterations.
m = lbfgsparams.m;  % Set memory parameter.
                    % NOTE:  Practical experience indicates that modest values
                    %        of m (e.g. between 3 and 20) often produce
                    %        satisfactory results.
xc.p = x.p;  % Set the current point to the initial point, x.p.

%  Initialize parameter structure for StepSizeSW function call.
params = struct('ftol', 1e-4, 'gtol', 0.9, 'xtol', 1e-6, 'stpmin', 0, ...
                'stpmax', 1e20, 'maxfev', 10000);

n = size(xc.p, 1);  % Dimension of vector space.

%  Allocate alpha and rho arrays for two-loop recursion procedure, and
%  allocate s and y arrays for storing a modified inverse Hessian approximation
%  implicitly via {s_i,y_i} pairs.  The rho, s, and y arrays will be used as
%  queues, with the m most recent entries being stored in the queues.
alfa_arr = zeros(1, m);
rho_arr = zeros(1, m);
s_arr = zeros(n, m);
y_arr = zeros(n, m);

I = eye(n);  % Locally stored identity matrix.

for i = 1:maxit
    %  Compute function and gradient at current point.
    xc.f = feval(func, xc.p, 1);
    xc.g = feval(func, xc.p, 2);
    
    %  Check for termination condition: (scaled) norm of gradient less
    %  than toler.
    if norm(xc.g) / min(1000, 1+abs(xc.f)) < toler
        inform.status = 1;  % Indicates success.
        inform.iter = i;  % Number of iterations.
        x.p = xc.p;
        x.f = xc.f;
        return;
    end
    
    %  For the first step, we use the identity matrix as an initial inverse
    %  Hessian approximation.
    if i == 1
        p = -xc.g;
    else
        %  Otherwise, we update the current inverse Hessian approximation, H,
        %  using the L-BFGS two-loop recursion procedure detailed in Algorithm
        %  7.4 of Nocedal and Wright.
        
        k = max(1, m-i+2);  % Sets (inclusive) lower index bound.  If
                            % i-1 >= m, we use m {s_j,y_j} vector pairs to
                            % obtain the approximate inverse Hessian and
                            % gradient product; otherwise, we use i-1 {s_j,y_j}
                            % vector pairs.
                            
        %  Scaling factor, gamma, attempts to estimate the size of the true
        %  Hessian matrix along the most recent search direction;
        %  calculated via (7.20) of Nocedal and Wright.
        gamma = (s_arr(:, k)'*y_arr(:,k)) / ((y_arr(:,k)'*y_arr(:,k)));
        %  The following choice for an initial Hessian approximation has proved
        %  effective in practice.
        H_init = gamma * I;
        q = xc.g;
        %  Compute product of inverse Hessian approximation and gradient for
        %  current iteration, to be stored in r.  See Algorithm 7.4 of
        %  Nocedal and Wright.
        for j = m : -1 : k
            alfa_arr(1,j) = rho_arr(1,j) * s_arr(:,j)' * q;
            q = q - alfa_arr(1,j) * y_arr(:,j);
        end
        r = H_init * q;
        for j = k : m
            beta = rho_arr(1, j) * y_arr(:,j)' * r;
            r = r + s_arr(:, j) * (alfa_arr(1, j) - beta);
        end;
        p = -r;  % Set search direction to the negative of this product.
    end
    
    %  Get step size that satisfies strong Wolfe conditions.
    %  NOTE:  alfa = 1 should always be tried first since this step length will
    %         eventually always be accepted (under certain conditions), thereby
    %         producing superlinear convergence of the overall algorithm.
    %         See page 142 of Nocedal and Wright.
    [alfa, x] = StepSizeSW(func, xc, p, 1, params);
    %  Update current point in p-direction with step size alpha.
    xc.p = xc.p + alfa * p;
    
    %  Update parameters.
    s = alfa * p;  %  s = new_point - prev_point = alfa * p
    y = feval(func, xc.p, 2) - xc.g;  % y = grad(new_point) - grad(prev_point)
    rho = 1 / (y'*s);
    
    %  Now we dequeue the oldest entry in each of the rho, s, and y queues, and
    %  enqueue the most recent rho, s, and y values.

    %  Move each rho, s, and y entry up one spot in their respective arrays.
    for j = 1 : m-1
        rho_arr(1, j) = rho_arr(1, j+1);
        s_arr(:, j) = s_arr(:, j+1);
        y_arr(:, j) = y_arr(:, j+1);
    end
    %  Set the the last entry in the rho, s, and y arrays to the most recent
    %  rho, s, and y values, respectively.
    rho_arr(1, m) = rho;
    s_arr(:, m) = s;
    y_arr(:, m) = y;
end
%  If reached, method failed.
inform.status = 0;  % Update status to failure indicator, 0.
inform.iter = maxit;  % Number of iterations.
x.p = xc.p;
x.f = feval(func, x.p, 1);
x.g = feval(func, x.p, 2);
return;  % Return inform and final point x
end