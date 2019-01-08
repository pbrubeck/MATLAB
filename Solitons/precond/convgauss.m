function [u] = convgauss(sig,B,jkrm,r,u)
u=hank(B,jkrm,u);
u=diag(exp(-(sig*r).^2/4))*u;
u=ihank(B,jkrm,u);
end