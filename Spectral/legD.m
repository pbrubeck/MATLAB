function [D, x, w] = legD(N)
% Legendre spectral differentiation matrix
k=1:N-1;
E=k./sqrt(4*k.*k-1);
[x,V]=trideigs(zeros(1,N), E); 
w=2*V(1,:).^2;
x=x(:); w=w(:);

X=repmat(x, [1, N]);
p=((-1).^(0:N-1)./abs(V(1,:)))'./sqrt(1-x.^2);
dX=X-X'+eye(N);
D=(p*(1./p)')./dX-diag(sum(1./dX));
end