function [inform, x] = GaussN(fun, resid, x, gnparams)
%  Implements the Gauss-Newton algorithm for solving non-linear least squares
%  problems.
%
%  Input:
%
%    fun      - a pointer to a function
%    resid    - a pointer to residual function
%    x        - the following structure:
%               * x.p - the starting point values
%               * x.f - the function value of x.p
%               * x.g - the gradient value of x.p
%               * x.h - the Hessian value of x.p
%               with only x.p guaranteed to be set.
%    gnparams - the following structure, as an example:
%
%                 gnparams = struct('maxit',1000,'toler',1.0e-4,
%                                  'lsmethod','chol');
%              
%              Note that 'lsmethod' can be 'chol', 'qr', or 'svd'.
%
%  Output:
%
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

%  Populate local caching of gnparams parameters.
toler = gnparams.toler;  % Set gradient tolerance.
maxit = gnparams.maxit;  % Set maximum number of allowed iterations.
lsmethod = gnparams.lsmethod;  % Set method to find search direction.

%  Initialize parameter structure for StepSizeSW function call.
params = struct('ftol', 1e-4, 'gtol', 0.9, 'xtol', 1e-6, 'stpmin', 0, ...
                'stpmax', 1e20, 'maxfev', 10000);

alfa = 1;  % Initial value of alpha.
xc.p = x.p;  % Set the current point to the initial point, x.p.

for i = 1:maxit
    %  Compute function and gradient at current point.
    xc.f = feval(fun, xc.p, 1);
    xc.g = feval(fun, xc.p, 2);
    
    %  Check for termination condition: norm of gradient less than toler.
    if norm(xc.g) < toler
        inform.status = 1;  % Indicates success.
        inform.iter = i;  % Number of iterations.
        x.p = xc.p;
        x.f = xc.f;
        x.g = xc.g;
        return;
    end

    %  Calculate residual function value and Jacobian at current point.
    r = feval(resid, xc.p, 1);
    J = feval(resid, xc.p, 2);

    %  Determine which method should be used to find the search direction, p.
    switch lsmethod
        case 'chol'  % Use Cholesky factorization of J'*J to get p.
            %  If J'*J is positive-definite, R upper triangular matrix such
            %  that R'*R = J'*J.
            %  Else, Cholesky failed and flag is a positive integer.
            [R, flag] = chol(J'*J);
            %  If the Cholesky factorization failed, return with fail status.
            if flag ~= 0
                inform.status = 0;  % Update status to failure indicator, 0.
                inform.iter = i;  % Number of iterations.
                x.p = xc.p;
                x.f = xc.f;
                return;
            end
            %  J'Jp = -J'r and R'R = -J'r, hence p = -R\(R'\(J'*r)).
            p = -R \ (R' \ (J'*r));
        case 'qr'  % Use QR factorization of J to get p.
            %  P permutation matrix, Q unitay matrix, and R upper triangular
            %  matrix with diagonal arranged in absolute decreasing order.
            [Q, R, P] = qr(J);
            n = size(J, 2);
            Q1 = Q(1:end, 1:n);
            R = R(1:n, 1:end);
            %  p = argmin ||J*p+r||^2 = solution of R*P'*p + Q1'*r = 0, hence
            %  p = -P' \ (R \ (Q1'*r)) = -P * (R \ (Q1' * r)).
            p = -P * (R \ (Q1'*r));
        case 'svd'  % Use SVD factorization of J to get p.
            %  U and V unitary matrices, S diagonal matrix.
            [U, S, V] = svd(full(J));
            n = size(J, 2);
            U1 = U(1:end, 1:n);
            S = S(1:n, 1:n);
            %  Since the Moore-Penrose inverse is pinv(J) = V*inv(S)*U1' and 
            %  p = -pinv(J)*r, we have that p = -V*inv(S)*U1'*r.
            p = -V * inv(S) * U1' * r;
    end
    
    %  Get step size that satisfies strong Wolfe conditions.
    [alfa, x] = StepSizeSW(fun, xc, p, alfa, params);
    %  Update current point in p-direction with step size alpha.
    xc.p = xc.p + alfa * p;
end
%  If reached, method failed.
inform.status = 0;  % Update status to failure indicator, 0.
inform.iter = maxit;  % Number of iterations i = maxit at this point.
x.p = xc.p;
x.f = feval(fun, x.p, 1);
x.g = feval(fun, x.p, 2);
return;  % Return inform and final point x
end
