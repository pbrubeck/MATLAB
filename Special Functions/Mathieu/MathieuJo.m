function [jom] = MathieuJo(m, q, z)
% Odd Radial Mathieu Function
n=42;
B=MathieuB(m, q, n);

sigma=(-1).^(0:n-1);
v1=sqrt(q)*exp(-z);
v2=sqrt(q)*exp(z);
if mod(m,2)==0
    s=((2:2:2*n)*B)*((sigma.*(2:2:2*n))*B);
    jom=-s/B(1)^2/q*bessprod3(B.*sigma', v1, v2);
else
    s=((1:2:2*n-1)*B)*(sigma*B);
    jom=s/B(1)^2/sqrt(q)*bessprod4(B.*sigma', v1, v2);
end
end

function [bp]=bessprod3(A, v1, v2)
a=besselj(0,v1);
b=besselj(0,v2);
aa=besselj(1,v1);
bb=besselj(1,v2);
bp=zeros(size(v1));
for j=1:length(A)
    aaa=2*(j)*aa./v1-a;
    bbb=2*(j)*bb./v2-b;
    bp=bp+A(j)*(aaa.*b-a.*bbb);
    [a, aa]=deal(aa, aaa);
    [b, bb]=deal(bb, bbb);
end
end

function [bp]=bessprod4(A, v1, v2)
a=besselj(0,v1);
b=besselj(0,v2);
aa=besselj(1,v1);
bb=besselj(1,v2);
bp=A(1)*(a.*bb-aa.*b);
for j=2:length(A)
    [a, aa]=deal(aa, 2*(j-1)*aa./v1-a);
    [b, bb]=deal(bb, 2*(j-1)*bb./v2-b);
    bp=bp+A(j)*(a.*bb-aa.*b);
end
end