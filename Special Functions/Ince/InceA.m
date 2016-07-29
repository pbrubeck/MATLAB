function [A] = InceA(p, m, q)
% Fourier coefficients for even Ince Polynomials (InceC)
if any((-1).^(p-m)==-1); error('p and m must have same parity'); end;
s=mod(p,2); k=(m-s)/2+1; n=(p-s)/2+1;
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
A=bsxfun(@times, A, sign(sum(A))./sqrt(1+(1-s)*A(1,:).^2));
A=A(:,order(k));
end