function [A] = InceA(p, m, q)
% Fourier coefficients for even Ince Polynomials (InceC)
n=floor(p/2)+1;
s=mod(m,2); k=(m-s)/2+1;
s=mod(p,2);
if s==0
    d=4*(0:n-1).^2;
    e1=q*(p/2+(1:n-1));
    e2=[q*p, q*(p/2-(1:n-2))];
else
    d=[1+q*(p+1)/2, (2*(1:n-1)+1).^2];
    e1=q*((p+1)/2+(1:n-1));
    e2=q*((p+1)/2-(1:n-1));
end

M=diag(d)+diag(e1,1)+diag(e2,-1);
[A,eta]=eig(M,'vector');
[eta,order]=sort(eta);
A=A(:,order(k));
A=A*sign(sum(A))/sqrt((1-s)*A(1)^2+1);
end