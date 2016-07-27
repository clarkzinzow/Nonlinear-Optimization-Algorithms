function [inform,x] = DogLeg(fun,x,nparams)
%  DogLeg implements the dogleg method for finding a solution to the subproblem
%
%       min m(p) = f + g'p + 1/2 * p' B p        s.t. ||p|| <= Del
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
%                           'maxdel',100,'eta',0.1,'method','chol',
%                           'hessian','exact','fail','cauchy');
%
%  Output:
%    inform - structure containing two fields:
%      * inform.status - 1 if gradient tolerance was
%                        achieved;
%                        0 if not.
%      * inform.iter   - the number of steps taken
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

for i = 1:maxit
    xc.f = feval(fun, xc.p, 1);  % Compute function at current point.
    xc.g = feval(fun, xc.p, 2);  % Compute gradient at current point.
    xc.h = sparse(feval(fun, xc.p, 4));  % Compute Hessian at current point.
    
    %  Check check for termination condition: norm of gradient less than toler.
    if norm(xc.g) < toler
        inform.status = 1;  % Indicates success.
        inform.iter = i;  % Number of iterations.
        x.p = xc.p;
        x.f = xc.f;
        x.g = xc.g;
        x.h = xc.h;
        return;
    end
    
    %  Calculate the Cauchy point, pU.
    
    %  If g'Hg <= 0 ...
    if xc.g'*xc.h*xc.g <= 0
        tau = 1;  % Set tao = 1.
    %  Otherwise (g'Hg > 0) ...
    else
        %  tau = min(||g||^3 / (del * g'Hg),1)
        tau = min(norm(xc.g)^3 / (del*xc.g'*xc.h*xc.g),1);
    end
    pU = -tau * del * xc.g / norm(xc.g);
    
    %  Calculate dogleg point, p.
    
    % Get the Cholesky factorization of the Hessian at the current point.
    [R,flag] = chol(xc.h);
    numFact = numFact + 1;  % Increment the 'chol' call counter.
    %  If the Cholesky factorization failed and nparams.fail = 'cauchy',
    %  then the dogleg point is the Cauchy point, pU.
    if flag ~= 0
        if isfield(nparams, 'fail') && strcmp(nparams.fail, 'cauchy')
            p = pU;  % Set dogleg point to Cauchy point.
        else
            %  Otherwise, repeatedly try 'multiple of the identity' Hessian
            %  modifications until the resulting matrix is positive-definite.
            beta = 0.001;  % t-shift parameter.
            minDiag = min(diag(xc.h));  % Minimum of the diagonal of Hessian.
            %  If diag > 0 componentwise...
            if minDiag > 0
                t = 0;  % Set initial shift to 0;
            else  % Otherwise..
                %  Set the initial shift to the negative of minDiag, plus the
                %  beta shift parameter (so diag > 0 after this shift is applied.
                t = -minDiag + beta;
            end
            %  Keep attempting Cholesky factorizations, with increasing shifts,
            %  until successful.
            while flag ~= 0
                [R,flag] = chol(xc.h + t*eye(size(xc.h,1)));
                numFact = numFact + 1;
                t = max(2*t,beta);
            end
            %  Get following solution.
            pB = -R \ (R' \ xc.g);
            %  Dogleg point is at the boundary of pU, pC-pU, and del.
            p = boundary(pU, pB - pU, del);
%           a = norm(pU)^2;
%           b = 2*pU'*(pB-pU);
%           c = norm(pU)^2 - del^2;
%           alpha = (-b + sqrt(b^2 - 4*a*c)) / (2*a);
%           p = pU + alpha * (pB-pU);
            %p = pU + ((-(2*pU'*(pB-pU)) + sqrt((2*pU'*(pB-pU))^2 - ...
            %4*norm(pU)^2*(norm(pU)^2-del^2))) / (2*norm(pU)^2)) * (pB - pU);
        end
    else
        %  Get following solution.
        pB = -R \ (R' \ xc.g);
        %  Dogleg point is at the boundary of pU, pC-pU, and del.
        p = boundary(pU, pB - pU, del);
        %a = norm(pU)^2;
        %b = 2*pU'*(pB-pU);
        %c = norm(pU)^2 - del^2;
        %alpha = (-b + sqrt(b^2 - 4*a*c)) / (2*a);
        %p = pU + alpha * (pB-pU);
        %p = pU + ((-(2*pU'*(pB-pU)) + sqrt((2*pU'*(pB-pU))^2 - ...
        %4*norm(pU)^2*(norm(pU)^2-del^2))) / (2*norm(pU)^2)) * (pB - pU);
    end
    
    %  Compute the reduction ratio.
    rho = (xc.f - feval(fun,xc.p + p,1)) / (-xc.g'*p - 0.5*p'*xc.h*p);
    
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
x.p = xc.p;
x.f = xc.f;
x.g = feval(fun, x.p, 2);
x.h = sparse(feval(fun,x.p,4));
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