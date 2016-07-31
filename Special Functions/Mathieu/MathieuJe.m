function [jem] = MathieuJe(m, q, z)
% Even Radial Mathieu Function Je_m(q; z)
if q==0
    jem=cosh(m*z);
    return;
end
n=42; s=mod(m,2);
A=MathieuA(m, q, n);
AA=A; AA(2:2:end)=-AA(2:2:end);
v1=sqrt(q)*exp(-z);
v2=sqrt(q)*exp(z);
if s==0
    p=sum(A)*sum(AA);
    jem=p/A(1)^2*bessprod1(AA, v1, v2);
else
    p=sum(A)*((1:2:2*n-1)*AA);
    jem=p/A(1)^2/sqrt(q)*bessprod2(AA, v1, v2);
end
end

function [bp]=bessprod1(A, v1, v2)
bp=zeros(size(v1));
j=0; tol=eps; r=1;
while(j<length(A) && any(r(:)>tol))
    delta=A(j+1)*(besselj(j,v1).*besselj(j,v2));
    bp=bp+delta;
    r=abs(delta)./(abs(bp)+tol);
    j=j+1;
end
end

function [bp]=bessprod2(A, v1, v2)
a0=besselj(0,v1);
b0=besselj(0,v2);
a1=besselj(1,v1);
b1=besselj(1,v2);
bp=A(1)*(a0.*b1+a1.*b0);
j=2; tol=1e-10; r=1;
while(j<=length(A) && any(r(:)>tol))
    [a0, a1]=deal(a1, besselj(j,v1));
    [b0, b1]=deal(b1, besselj(j,v2));
    delta=A(j)*(a0.*b1+a1.*b0);
    bp=bp+delta;
    r=abs(delta)./(abs(bp)+tol);
    j=j+1;
end
end