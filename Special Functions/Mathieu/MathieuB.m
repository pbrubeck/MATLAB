function [B] = MathieuB(m, q, n)
% Fourier coefficients for odd Mathieu Functions (MathieuS, MathieuJo)
E(1:n-1)=q;
if mod(m,2)==0
    D=(2:2:2*n).^2;
else
    D=[1-q, (3:2:2*n-1).^2];  
end
[~,B]=trideigs(D,E);
B=B(:,ceil(m/2));
end

