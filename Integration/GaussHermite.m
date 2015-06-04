function [x, w] = GaussHermite(n)
    H=Hermite(n+1);
    x=polyRoots(H(n+1,:));
    w=(2^n*factorial(n+1)*sqrt(pi))./(n^2*(Horner(H(n,:), x)).^2);
end
