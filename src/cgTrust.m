function [inform, x] = cgTrust(fun, x, cgtparams)
%  Implements the Steihaug-Toint conjugate gradient trust region method for
%  finding an approximate solution to the subproblem:
%
%        min m(p) = f + g'p + 1/2 * p' B p        s.t. ||p|| <= Del
%
%  Input:
%    fun      - a pointer to a function
%    x        - the following structure:
%               * x.p - the starting point values
%               * x.f - the function value of x.p
%               * x.g - the gradient value of x.p
%               * x.h - the Hessian value of x.p
%               with only x.p guaranteed to be set.
%    cgtparams - the following structure, as an example:
%         cgtparams = struct('maxit',1000,'toler',1.0e-4,'initdel',1,
%                          'maxdel',100,'eta',0.1);
%
%  Output:
%    inform - structure containing two fields:
%      * inform.status - 1 if gradient tolerance was
%                        achieved;
%                        0 if not.
%      * inform.iter   - the number of steps taken
%      * inform.cgiter - the number of conjugate iterations
%    x      - the solution structure, with the solution point along with
%             function, gradient, and Hessian evaluations thereof.

%  Number of function, gradient, and Hessian evaluations, and number of Cholesky
%  factorizations.
global numf numg numh numFact
numf = 0;
numg = 0;
numh = 0;
numFact = 0;

%  Populate local caching of cgtparams parameters.
toler = cgtparams.toler;  % Set gradient tolerance.
maxit = cgtparams.maxit;  % Set maximum number of allowed iterations.
initdel = cgtparams.initdel;  % Set initial delta value.
maxdel = cgtparams.maxdel;  % Set maximum delta value.
eta = cgtparams.eta;  % Set eta value.

del = initdel;  % Set delta value to initial delta value.
xc.p = x.p;  % Set the current point to the initial point, x.p.
cgiter = 0;  % Number of conjugate iterations.

for i = 1:maxit
    %  Compute function, gradient, and Hessian at current point.
    xc.f = feval(fun, xc.p, 1);
    xc.g = feval(fun, xc.p, 2);
    xc.h = sparse(feval(fun, xc.p, 4));
    
    %  Check check for termination condition: norm of gradient less than toler.
    if norm(xc.g) < toler
        inform.status = 1;  % Indicates success.
        inform.iter = i;  % Number of iterations.
        inform.cgiter = cgiter;  % Number of conjugate iterations.
        x.p = xc.p;
        x.f = xc.f;
        x.g = xc.g;
        x.h = xc.h;
        return;
    end
    
    %  Conjugate gradient method.
    
    g = xc.g;  %  Set residual to gradient.
    initGradToler = toler * norm(g);  % Set residual norm-based toler.
    d = -g;  % Current search direction.
    k = 0;  % Iteration index.
    z = zeros(size(xc.h, 1), 1);
    while 1
        %  If current search direction, d, is a direction of nonpositive
        %  curvature along xc.h...
        if d' * xc.h * d <= 0
            %  Find the boundary point and break.
            p = boundary(z, d, del);
            break;
        end
        alpha = norm(g)^2 / (d' * xc.h * d);
        %  If z+alpha*d violates the trust-region bound...
        if norm(z + alpha * d) >= del
            %  Find the boundary point and break.
            p = boundary(z, d, del);
            break;
        end
        z = z + alpha * d;
        gN = g + alpha * xc.h * d;
        %  If current gradient iterate is below threshold...
        if or(norm(gN) < toler, norm(gN) < initGradToler)
            %  Set point to current z iterate and break.
            p = z;
            break;
        end
        
        %  Parameter updates.
        beta = norm(gN)^2 / norm(g)^2;  % Update beta.
        d = -gN + beta * d;  % Update search direction.
        g = gN;  % Update gradient.
        k = k + 1;  % Update iteration index.    
    end
    cgiter = cgiter + k;  % Update number of conjugate gradient iterations.
    
    %  Compute the reduction ratio.
    rho = (xc.f - feval(fun,xc.p + p,1))/(-xc.g'*p - 0.5*p'*xc.h*p);
    
    %  Update the trust region; i.e., update del and the current point.
    if rho < 0.25
        del = del / 4;
    else
        %  Note that radius only increases if norm of p reaches the trust
        %  region boundary.
        if rho > 0.75 && norm(p) == del
            del = min(2*del, maxdel);
        end
    end
    if rho > eta
        xc.p = xc.p + p;
    end
end
%  If reached, method failed.
inform.status = 0;  % Update status to failure indicator, 0.
inform.iter = maxit;  % Number of iterations = i = maxit at this point.
inform.cgiter = cgiter;  % Number of conjugate iterations.
x.p = xc.p;
x.f = feval(fun, x.p, 1);
x.g = feval(fun, x.p, 2);
x.h = sparse(feval(fun, x.p, 4));
return;  % Return inform and final point x
end

function p = boundary(p, q, del)
%  Finds the boundary point.
a = norm(q)^2;
b = 2*p'*q;
c = norm(p)^2 - del^2;
alpha = (-b + sqrt(b^2 - 4*a*c)) / (2*a);
p = p + alpha * q;
return;
end