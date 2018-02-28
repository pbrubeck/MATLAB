function [C] = legC(xleg, x)
% Cardinal Legendre interpolation matrix
N=length(xleg);
xleg=xleg(:);
x=x(:);
c=[zeros(N-1,1);1];
P0=LegendreP(c,0,xleg);
P1=LegendreP(c,1,x);
C0=repmat(sqrt(1-x.^2).*P1, [1,N]);
C1=repmat(N*(N-1)*P0', [length(x),1]);
X=repmat(x, [1,N]);
Xj=repmat(xleg', [length(x),1]);
dX=X-Xj;
C=(C0./C1)./dX;
C(abs(dX)<10*eps)=1;
end