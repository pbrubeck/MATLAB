function [C] = legC(xleg, x)
% Cardinal Legendre interpolation matrix
N=length(xleg);
xleg=xleg(:);
x=x(:);
c=[zeros(N-1,1);1];
C0=sqrt(1-x.^2).*LegendreP(c,1,x);
C1=N*(N-1)*LegendreP(c,0,xleg);
X=repmat(x, [1,N]);
Xj=repmat(xleg', [length(x),1]);
dX=X-Xj;
C=(C0*(1./C1)')./dX;
C(abs(dX)<10*eps)=1;
%C=diag(1./sum(C,2))*C;
end