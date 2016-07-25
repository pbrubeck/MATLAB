function [jem] = MathieuJe(m, q, z)
% Even Radial Mathieu Function
if q==0
    jem=cosh(m*z);
    return;
end
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
a0=besselj(0,v1);
b0=besselj(0,v2);
a1=besselj(1,v1);
b1=besselj(1,v2);
bp=A(1)*(a0.*b0);
for j=2:length(A)
    [a0, a1]=deal(a1, besselj(j,v1)); %2*(j-1)*a1./v1-a0);
    [b0, b1]=deal(b1, besselj(j,v2)); %2*(j-1)*b1./v2-b0);
    bp=bp+A(j)*(a0.*b0);
end
end

function [bp]=bessprod2(A, v1, v2)
a0=besselj(0,v1);
b0=besselj(0,v2);
a1=besselj(1,v1);
b1=besselj(1,v2);
bp=A(1)*(a0.*b1+a1.*b0);
for j=2:length(A)
    [a0, a1]=deal(a1, besselj(j,v1)); %2*(j-1)*a1./v1-a0);
    [b0, b1]=deal(b1, besselj(j,v2)); %2*(j-1)*b1./v2-b0);
    bp=bp+A(j)*(a0.*b1+a1.*b0);
end
end