function [inform, x] = GaussN(fun, x, nparams)
% GaussN implements the Gauss Newton method for finding a solution to
%
%       min f(x) = 1/2 * sum_{j=1}^m r_j^2(x)
%
%  Input:
%
%    fun      - a pointer to a function
%    x        - the following structure:
%               * x.p - the starting point values
%               * x.f - the function value of x.p
%               * x.g - the gradient value of x.p
%               * x.h - the Hessian value of x.p
%               with only x.p guaranteed to be set.
%    nparams - the following structure, as an example:
%
%                 nparams = struct('maxit',1000,'toler',1.0e-4,
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

global numf numg numh resid  % resid - function pointer defined in hwk6.m.

numf = 0;  % Number of function evaluations.
numg = 0;  % Number of gradient evaluations.
numh = 0;  % Number of Hessian evaluations.

%  Populate local versions of nparams parameters.
toler = nparams.toler;  % Set gradient tolerance.
maxit = nparams.maxit;  % Set maximum number of allowed iterations.
lsmethod = nparams.lsmethod;  % Set method to find search direction.

%  Populate params structure, to be provided to StepSize.
params = struct('ftol', 1e-4, 'gtol', 0.9, 'xtol', 1e-6, 'stpmin', 0, ...
                'stpmax', 1e20, 'maxfev', 10000);

alfa = 1;  % Initial value of alpha.
xc.p = x.p;  % Set the current point to the initial point, x.p.

for i = 1:maxit
    xc.f = feval(fun, xc.p, 1);  % Compute function at current point.
    xc.g = feval(fun, xc.p, 2);  % Compute gradient at current point.
    
    %  Check if norm of gradient at current point is below toler.
    if norm(xc.g) < toler
        inform.status = 1;  % Indicates success.
        inform.iter = i;  % Set number of iterations.
        x.p = xc.p;
        x.f = xc.f;
        return;
    end

    %  Calculate residual function value and Jacobian at current point.
    r = feval(resid, xc.p, 1);  % Residual function value at current point.
    J = feval(resid, xc.p, 2);  % Residual Jacobian at current point.

    %  Determine which method should be used for finding the search direction.
    switch lsmethod
        %  Get the Cholesky factorization of J'*J, use the result to get pgn.
        case 'chol'
            %  If J'*J is pd, R upper triangular matrix such that R'*R = J'*J.
            %  Else, Cholesky failed and flag is a positive integer.
            [R, flag] = chol(J'*J);
            %  If the Cholesky factorization failed, return with fail status.
            if flag ~= 0
                inform.status = 0;
                inform.iter = maxit;
                x.p = xc.p;
                x.f = xc.f;
                return;
            end
            %  J'Jp = -J'r and R'R = -J'r, hence pgn = -R\(R'\(J'*r)).
            pgn = -R \ (R' \ (J' * r));
        %  Get the QR factorization of J, use the result to get pgn.
        case 'qr'
            %  P permutation matrix, Q unitay matrix, and R upper triangular
            %  matrix with diagonal arranged in absolute decreasing order.
            [Q R P] = qr(J);
            n = size(J, 2);
            Q1 = Q(1:end, 1:n);
	    R = R(1:n, 1:end);
            %  By Q being unitary and P being a permutation matrix.
            pgn = -P * (R \ (Q1' * r));
        %  Get the SVD factorization of J, use the result to get pgn.
        case 'svd'
            %  U and V unitary matrices, S diagonal matrix.
            [U S V] = svd(full(J));
            n = size(J, 2);
            U1 = U(1:end, 1:n);
            S = S(1:n, 1:n);
            %  By results in lecture 34.
            pgn = -V * inv(S) * U1' * r;
    end
    
    %  Calculate the step size and the current candidate solution point.
    [alfa, x] = StepSizeSW(fun, xc, pgn, alfa, params);
    %  Update our current point by taking a step in the calculated direction.
    xc.p = xc.p + alfa * pgn;
    
end
%  If reached, method failed.
inform.status = 0;  % Update status to failure indicator, 0.
inform.iter = maxit;  % Number of iterations i = maxit at this point.
x.p = xc.p;
x.f = xc.f;
return;  % Return inform and final point x
end
