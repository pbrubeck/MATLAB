function [C] = chebC(xcheb, x)
% Cardinal Chebyshev interpolation matrix
N=length(xcheb);
xcheb=xcheb(:);
x=x(:);
P0=ChebT([zeros(N-1,1);1],xcheb);
P0([1,end])=2*P0([1,end]);
P1=ChebU([zeros(N-2,1);1],x);
C0=repmat((x.^2-1).*P1, [1,N]);
C1=repmat((N-1)*P0', [length(x),1]);
X=repmat(x, [1,N]);
Xj=repmat(xcheb', [length(x),1]);
dX=X-Xj;
C=(C0./C1)./dX;
C(abs(dX)<10*eps)=1;
end