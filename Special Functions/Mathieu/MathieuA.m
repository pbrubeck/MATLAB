function [A] = MathieuA(m, q, n)
% Fourier coefficients for even Mathieu Functions (MathieuC, MathieuJe)
E(1:n-1)=q;
s=mod(m,2); k=(m-s)/2+1;
if s==0
    D=(0:2:2*n-2).^2;
    E(1)=sqrt(2)*q;
else
    D=[1+q, (3:2:2*n-1).^2];  
end
[~,A]=trideigs(D,E);
A=bsxfun(@times, A, sign(sum(A)));
A(1,:)=A(1,:)*(1-(1-s)/(2+sqrt(2)));
A=A(:,k);
end