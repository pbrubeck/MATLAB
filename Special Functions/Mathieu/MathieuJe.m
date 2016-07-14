function [jem] = MathieuJe(m, q, z)
% Even Radial Mathieu Function
n=42;
A=MathieuA(m, q, n);

sigma=(-1).^(0:n-1);
v1=sqrt(q)*exp(-z);
v2=sqrt(q)*exp(z);
if mod(m,2)==0
    p=sum(A)*(sigma*A);
    jem=p/A(1)^2*bessprod1(A.*sigma', v1, v2);
else
    p=sum(A)*(sigma.*(1:2:2*n-1))*A;
    jem=p/A(1)^2/sqrt(q)*bessprod2(A.*sigma', v1, v2);
end
end

function [bp]=bessprod1(A, v1, v2)
a=besselj(0,v1);
b=besselj(0,v2);
aa=besselj(1,v1);
bb=besselj(1,v2);
bp=A(1)*(a.*b);
for j=2:length(A)
    [a, aa]=deal(aa, 2*(j-1)*aa./v1-a);
    [b, bb]=deal(bb, 2*(j-1)*bb./v2-b);
    bp=bp+A(j)*(a.*b);
end
end

function [bp]=bessprod2(A, v1, v2)
a=besselj(0,v1);
b=besselj(0,v2);
aa=besselj(1,v1);
bb=besselj(1,v2);
bp=A(1)*(a.*bb+aa.*b);
for j=2:length(A)
    [a, aa]=deal(aa, 2*(j-1)*aa./v1-a);
    [b, bb]=deal(bb, 2*(j-1)*bb./v2-b);
    bp=bp+A(j)*(a.*bb+aa.*b);
end
end