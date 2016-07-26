function [B] = InceB(p, m, q)
% Fourier coefficients for odd Ince Polynomials (InceS)
N=floor(p/2)+1;
s=2-mod(m,2); k=(m-s)/2+1;
if s==1
    d=[1-q*(p+1)/2, (2*(1:N-1)+1).^2];
    e1=q*((p+1)/2+(1:N-1));
    e2=q*((p+1)/2-(1:N-1));
else
    d=4*(1:N).^2;
    e1=q*(p/2+1+(1:N-1));
    e2=q*(p/2+1-(2:N));
end
    
M=diag(d)+diag(e1,1)+diag(e2,-1);
[B,eta]=eig(M,'vector');
[eta,order]=sort(eta);
B=B(:,order(k));
B=B*sign((1+s:2:2*N+s)*B);
end