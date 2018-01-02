function [D, x, w] = chebD(N)
% Computes the N by N Chebyshev differentiation matrix and grid.
x=cos(pi*(0:N-1)/(N-1))';
c=[2; ones(N-2,1); 2].*(-1).^(0:N-1)';
X=repmat(x,1,N);
dX=X-X';
D=(c*(1./c)')./(dX+eye(N));
D=D-diag(sum(D, 2));
end