function [lam,mu,x1,x2] = newton_mep(A1,B1,C1,A2,B2,C2,x1,x2,lam,mu)
% Multi-Parameter EigenProblem solved through Newton Raphson
% (A1-lam*B1-mu*C1)*x1=0
% (A2-lam*B2-mu*C2)*x2=0
m=size(A1,1);
n=size(A2,1);
J=zeros(m+n+2);
x=[normc(x1(:)); normc(x2(:)); lam; mu];
dx=ones(size(x));
while(max(abs(dx))>1E-10)
    x1=x(1:m);
    x2=x(m+1:end-2);
    lam=x(end-1);
    mu=x(end);
    
    D1=A1-lam*B1-mu*C1;
    D2=A2-lam*B2-mu*C2;
    F=[D1*x1; D2*x2; x1'*x1-1; x2'*x2-1];
    
    J(1:m,1:m)=D1;
    J(m+1:m+n,m+1:m+n)=D2;
    J(1:m+n,end-1:end)=[-B1*x1, -C1*x1; -B2*x2, -C2*x2];
    J(end-1,1:m)=2*x1';
    J(end,m+1:m+n)=2*x2';
    
    dx=J\F;
    x=x-dx;
end
end
