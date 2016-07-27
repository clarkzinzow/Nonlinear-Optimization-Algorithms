function rhoout = Rho(X,Y,alpha,beta,mode)
if bitand(mode,1)
    N = length(X);
    rhoout = zeros(N,1);
    for i=1:N
        rhoout(i) = 1+alpha*exp(-1*beta*(X(i).^2 + Y(i).^2));
    end
end
    return;
end