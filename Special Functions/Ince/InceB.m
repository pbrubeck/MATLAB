function [B] = InceB(p, m, q)
% Fourier coefficients for odd Ince Polynomials (InceS)
if any((-1).^(p-m)==-1); error('p and m must have same parity'); end;
s=2-mod(p,2); k=(m-s)/2+1; n=(p-s)/2+1;
if s==2
    d=4*(1:n).^2;
    e1=q*(p/2+1+(1:n-1));
    e2=q*(p/2+1-(2:n));
else
    d=[1-q*(p+1)/2, (2*(1:n-1)+1).^2];
    e1=q*((p+1)/2+(1:n-1));
    e2=q*((p+1)/2-(1:n-1));
end

M=diag(d)+diag(e1,1)+diag(e2,-1);
[B,eta]=eig(M,'vector');
[eta,order]=sort(eta);
B=bsxfun(@times, B, sign((1+s:2:2*n+s)*B));
B=B(:,order(k));
end