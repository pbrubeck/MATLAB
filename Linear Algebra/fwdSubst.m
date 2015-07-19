function x=fwdSubst(L, b)
% Solves Lx=b for lower triangular matrices.
n=size(L,1);
x=zeros(n,1);
for i=1:n
    x(i)=(b(i)-L(i,1:i-1)*x(1:i-1))/L(i,i);
end
end