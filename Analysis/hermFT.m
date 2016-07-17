function [F] = hermFT(x,w)
% Computes Fourier Transform operator given the Gauss-Hermite nodes and weights
n=size(x,1);
H=zeros(n);
a=1;
for i=1:n
    H(:,i)=HermitePsi(a,x(:));
    a=[0, a];
end
F=H*diag((1i).^(0:n-1))*H'*diag(w(:).*exp(x(:).^2));

F=1/sqrt(2*pi)*exp(-1i*(x*x'))*diag(w(:).*exp(x(:).^2));
end