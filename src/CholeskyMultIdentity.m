function [R, tau] = CholeskyMultIdentity(H)
%  Implements Cholesky with added multiple of the identity.  This attempts to
%  to find a scalar tau > 0 such that H + tau * I is sufficiently positive
%  definite, where I is the identity matrix.

global numFact  % Number of Cholesky factorizations attempted.

%  First, we intially try a Cholesky factorization.
[R, fail] = chol(H);
numFact = numFact + 1;
%  If it does not fail, we're done.
if fail == 0
    tau = 0;
    return;
end

%  If the initial Cholesky factorization fails, we attempt to find a scalar
%  tau > 0 such that H + tau * I is sufficiently positive definite.
beta = 0.001;  % Heuristic for increasing tau.
min_H_diag = min(diag(H));  % Smallest diagonal of H.

%  If smallest diagonal of H is positive, set tau to 0; otherwise, set to
%  nonnegative version of the smallest diagonal plus the beta heuristic.
if min_H_diag > 0
    tau = 0;
else
    tau = -min_H_diag + beta;
end

I = eye(size(H,1));  % Identity matrix.

%  Repeatedly add a tau-multiple of the identity to H until the Cholesky
%  factorization succeeds.  Upon each failure, double tau.
while 1
    [R, fail] = chol(H + tau * I);
    numFact = numFact + 1;
    if fail == 0
        return;
    else
        %  NOTE:  In order to decrease number of factorizations, we may want to
        %         increase tau by a factor of 10 instead of 2.
        tau = max(2*tau, beta);
    end
end
end
        