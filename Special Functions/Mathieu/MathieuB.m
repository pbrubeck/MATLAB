function [B] = MathieuB(m, q, n)
% Fourier coefficients for odd Mathieu Functions (MathieuS, MathieuJo)
E(1:n-1)=q;
s=2-mod(m,2); k=(m-s)/2+1;
if s==2
    D=(2:2:2*n).^2;
else
    D=[1-q, (3:2:2*n-1).^2];
end
[~,B]=trideigs(D,E);
B=bsxfun(@times, B, sign((1+s:2:2*n+s)*B));
B=B(:,k);
end