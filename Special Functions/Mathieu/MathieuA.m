function [A] = MathieuA(m, q, n)
% Fourier coefficients for even Mathieu Functions (MathieuC, MathieuJe)
E(1:n-1)=q;
if mod(m,2)==0
    D=(0:2:2*n-2).^2;
    E(1)=sqrt(2)*q;
else
    D=[1+q, (3:2:2*n-1).^2];  
end
[~,A]=trideigs(D,E);

k=1+floor(m/2);
A=A(:,k);
A=sign(A(k))*A;
A(1)=A(1)*(1-(1-mod(m,2))/(2+sqrt(2)));
end