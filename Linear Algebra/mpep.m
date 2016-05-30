function [u,v,lam,mu] = mpep(A,B,C,D,u,v,lam,mu)
% Multi-Parameter EigenProblem solved through Newton Raphson
% (A-mu*I)*u=lam*C*u, (B+mu*I)*v=lam*D*v
m=size(A,1);
n=size(B,1);
J=zeros(m+n+2);
x=[normc(u(:)); lam; mu; normc(v(:))];
dx=ones(size(x));
while(max(abs(dx))>1E-10)
    u=x(1:m);
    lam=x(m+1);
    mu=x(m+2);
    v=x(m+3:end);
    
    F=[u'*u-1; (A-mu*eye(m)-lam*C)*u; (B+mu*eye(n)-lam*D)*v; v'*v-1];
    J(1,1:m)=2*u';
    J(2:m+1,1:m+2)=[A-mu*eye(m)-lam*C, -C*u, -u];
    J(m+2:end-1,m+1:end)=[-D*v, v, B+mu*eye(n)-lam*D];
    J(end,m+3:end)=2*v';
    
    dx=J\F;
    x=x-dx;
end
end
