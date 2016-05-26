function [jem] = MathieuJe(m, n, q, z)
j=mod(m,2);
k=(m-j+2)/2;
E(1:n-1)=q;
if j==0
    D=(0:2:2*n-2).^2;
    E(1)=sqrt(2)*q;
else
    D=[1+q, (3:2:2*n-1).^2];  
end
[~,A]=trideigs(D,E);
A=A(:,k);
A=sign(A(k))*A;
A(1)=A(1)*(1-(1-j)/(2+sqrt(2)));
A=A.*(-1).^(0:n-1)';

v1=sqrt(q)*exp(-z);
v2=sqrt(q)*exp(z);
jem=bessprod(A, v1, v2)/(A(1)^2);
end

function [bp]=bessprod(A, v1, v2)
a=besselj(0,v1);
b=besselj(0,v2);
aa=besselj(1,v1);
bb=besselj(1,v2);
bp=A(1)*(a.*b);
for j=2:length(A)
    bp=bp+A(j)*(aa.*bb);
    [a,aa]=deal(aa, 2*(j-1)*aa./v1-a);
    [b,bb]=deal(bb, 2*(j-1)*bb./v2-b);
end
end