function [jom] = MathieuJo(m, q, z)
% Odd Radial Mathieu Function
if q==0
    jom=sinh(m*z);
    return;
end
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
a0=besselj(0,v1);
b0=besselj(0,v2);
a1=besselj(1,v1);
b1=besselj(1,v2);
bp=zeros(size(v1));
for j=1:length(A)
    a2=besselj(j+1,v1); %2*j*a1./v1-a0;
    b2=besselj(j+1,v2); %2*j*b1./v2-b0;
    bp=bp+A(j)*(a2.*b0-a0.*b2);
    [a0, a1]=deal(a1, a2);
    [b0, b1]=deal(b1, b2);
end
end

function [bp]=bessprod4(A, v1, v2)
a0=besselj(0,v1);
b0=besselj(0,v2);
a1=besselj(1,v1);
b1=besselj(1,v2);
bp=A(1)*(a0.*b1-a1.*b0);
for j=2:length(A)
    [a0, a1]=deal(a1, besselj(j,v1)); %2*(j-1)*a1./v1-a0);
    [b0, b1]=deal(b1, besselj(j,v2)); %2*(j-1)*b1./v2-b0);
    bp=bp+A(j)*(a0.*b1-a1.*b0);
end
end