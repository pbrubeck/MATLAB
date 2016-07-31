function [jom] = MathieuJo(m, q, z)
% Odd Radial Mathieu Function Jo_m(q; z)
if q==0
    jom=sinh(m*z);
    return;
end
n=42; s=2-mod(m,2);
B=MathieuB(m, q, n);
BB=B; BB(2:2:end)=-BB(2:2:end);
v1=sqrt(q)*exp(-z);
v2=sqrt(q)*exp(z);
if s==2
    p=((2:2:2*n)*B)*((2:2:2*n)*BB);
    jom=-p/B(1)^2/q*bessprod3(BB, v1, v2);
else
    p=((1:2:2*n-1)*B)*sum(BB);
    jom=p/B(1)^2/sqrt(q)*bessprod4(BB, v1, v2);
end
end

function [bp]=bessprod3(A, v1, v2)
a0=besselj(0,v1);
b0=besselj(0,v2);
a1=besselj(1,v1);
b1=besselj(1,v2);
bp=zeros(size(v1));
j=1; tol=1e-10; r=1;
while(j<=length(A) && any(r(:)>tol))
    a2=besselj(j+1,v1);
    b2=besselj(j+1,v2);
    delta=A(j)*(a2.*b0-a0.*b2);
    bp=bp+delta;
    r=abs(delta)./(abs(bp)+tol);
    [a0, a1]=deal(a1, a2);
    [b0, b1]=deal(b1, b2);
    j=j+1;
end
end

function [bp]=bessprod4(A, v1, v2)
a0=besselj(0,v1);
b0=besselj(0,v2);
a1=besselj(1,v1);
b1=besselj(1,v2);
bp=A(1)*(a0.*b1-a1.*b0);
j=2; tol=1e-10; r=1;
while(j<=length(A) && any(r(:)>tol))
    [a0, a1]=deal(a1, besselj(j,v1));
    [b0, b1]=deal(b1, besselj(j,v2));
    delta=A(j)*(a0.*b1-a1.*b0);
    bp=bp+delta;
    r=abs(delta)./(abs(bp)+tol);
    j=j+1;
end
end