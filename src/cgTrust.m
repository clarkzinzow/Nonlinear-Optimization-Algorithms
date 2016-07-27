function [inform, x] = cgTrust(fun, x, nparams)
%  cgTrust implements the Steihaug-Toint conjugate gradient trust region method
%  for finding an approximate solution to the subproblem:
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
%    nparams - the following structure, as an example:
%         nparams = struct('maxit',1000,'toler',1.0e-4,'initdel',1,
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

global numf numg numh numFact

numf = 0;  % Number of function evaluations.
numg = 0;  % Number of gradient evaluations.
numh = 0;  % Number of Hessian evaluations.
numFact = 0;  % Number of Cholesky factorizations.

%  Populate local versions of nparams parameters.
toler = nparams.toler;  % Set gradient tolerance.
maxit = nparams.maxit;  % Set maximum number of allowed iterations.
initdel = nparams.initdel;  % Set initial delta value.
maxdel = nparams.maxdel;  % Set maximum delta value.
eta = nparams.eta;  % Set eta.

del = initdel;  % Set delta value to initial delta value.
xc.p = x.p;  % Set the current point to the initial point, x.p.
cgiter = 0;  % Number of conjugate iterations.

for i = 1:maxit
    xc.f = feval(fun, xc.p, 1);  % Compute function at current point.
    xc.g = feval(fun, xc.p, 2);  % Compute gradient at current point.
    xc.h = sparse(feval(fun, xc.p, 4));  % Compute Hessian at current point.
    
    %  Check check for termination condition: norm of gradient less than toler.
    if norm(xc.g) < toler
        inform.status = 1;  % Indicates success.
        inform.iter = i;  % Set number of iterations.
        inform.cgiter = cgiter;  % Set number of conjugate iterations.
        x.p = xc.p;
        x.f = xc.f;
        x.g = xc.g;
        x.h = xc.h;
        return;
    end
    
    %  Conjugate gradient method.
    
    g = xc.g;  %  Set residual to gradient.
    initGradToler = toler * norm(g);  % Set residual norm-based toler.
    d = -g;
    k = 0;
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
        %  Update beta.
        beta = norm(gN)^2 / norm(g)^2;
        %  Update search direction.
        d = -gN + beta * d;
        %  Update gradient.
        g = gN;
        %  Update iteration index.
        k = k + 1;      
    end
    cgiter = cgiter + k;
    
    %  Compute the reduction ratio.
    rho = (xc.f - feval(fun,xc.p + p,1))/(-xc.g'*p - 0.5*p'*xc.h*p);
    
    %  Update the trust region; i.e., update del and the current point.

    if rho < 0.25
        del = del/4;
    else
        %  Note that radius only increases if ||p|| reaches the TR boundary.
        if rho > 0.75 && norm(p) == del
            del = min(2*del, maxdel);
        end
        % else del_{k+1} = del_k; i.e., del does not change.
    end
    if rho > eta
        xc.p = xc.p + p;
    end
    % else p_{k+1} = p_k; i.e., the current point does not change.
end
%  If reached, method failed.
inform.status = 0;  % Update status to failure indicator, 0.
inform.iter = maxit;  % Number of iterations i = maxit at this point.
inform.cgiter = cgiter;  % Number of conjugate iterations.
x.p = xc.p;
x.f = xc.f;
x.g = feval(fun, x.p, 2);
x.h = sparse(feval(fun,x.p,4));
return;  % Return inform and final point x
end

function p = boundary(p, q, del)
%  boundary function finds the boundary point.
a = norm(q)^2;
b = 2*p'*q;
c = norm(p)^2-del^2;
alpha = (-b+sqrt(b^2-4*a*c)) / (2*a);
p = p+alpha * q;
return;
end