function [] = MathieuTest(m, q)
L=25;
z=linspace(-L,L,10000);

n=100;
A=MathieuA(m, q, n);

sigma=(-1).^(0:n-1);
v1=sqrt(q)*exp(-z);
v2=sqrt(q)*exp(z);
s=mod(m,2);
if s==0
    AA=sigma(:).*A;
    p=sum(A)*sum(AA);
    [jem0, j]=bpexact(AA, v1, v2);
    jem0=p/A(1)^2*jem0;
    jem1=p/A(1)^2*bpclen(AA, v1, v2);
    jem2=p/A(1)^2*bprecursion(AA(1:j), v1, v2);
    jem3=p/A(1)^2*bpfactor(AA(1:j), v1, v2);
    jem4=sum(A)/A(1)*sumj2k(A(1:j), s, v2-v1);
    jem5=sum(AA)/A(1)*sumj2k(AA(1:j), s, v2+v1);
    
    figure(1);
    plot(z, [jem0-jem5]); title(sprintf('%d terms', j));
    threeTerm(v1, v2);
end


end

function []=threeTerm(q, z)
z=z(:);
y0=besselj(0,z);
y1=besselj(1,z);
r1=zeros(2, length(y0)-1);
r2=zeros(2, length(y0)-1);
for j=2:9
    y2=besselj(j,z);
    for i=1:length(y0)-1
        r1(:,i)=[y1(i:i+1), y0(i:i+1)]\y2(i:i+1);
    end
    [y0, y1]=deal(y1, y2);
    figure(2); plot(r1'); drawnow; pause(1);
end

end

function [bp, j]=bpexact(A, v1, v2)
bp=zeros(size(v1));
j=0; r=1; tol=eps;
while(j<=length(A) && any(r(:)>tol))
    delta=A(j+1)*(besselj(j,v1).*besselj(j,v2));
    bp=bp+delta;
    r=abs(delta)./(abs(bp)+tol);
    j=j+1;
end
end

function [bp]=bprecursion(A, v1, v2)
a0=besselj(0,v1);
b0=besselj(0,v2);
a1=besselj(1,v1);
b1=besselj(1,v2);
bp=A(1)*(a0.*b0);
for j=2:length(A)
    [a0, a1]=deal(a1, 2*(j-1)*a1./v1-a0);
    [b0, b1]=deal(b1, 2*(j-1)*b1./v2-b0);
    bp=bp+A(j)*(a0.*b0);
end
end

function [bp]=bpfactor(A, v1, v2)
bp=zeros(size(v1));
p0=zeros(size(v1));
p1=zeros(size(v1));
q0=zeros(size(v1));
q1=zeros(size(v1));
for alpha=0:1
    for beta=0:1
        p0(1:end)=1-alpha;
        p1(1:end)=alpha;
        q0(1:end)=1-beta;
        q1(1:end)=beta;
        S=A(1)*(p0.*q0);
        for j=2:length(A)
            [p0, p1]=deal(p1, 2*(j-1)*p1./v1-p0);
            [q0, q1]=deal(q1, 2*(j-1)*q1./v2-q0);
            S=S+A(j)*(p0.*q0);
        end
        bp=bp+S.*besselj(alpha, v1).*besselj(beta, v2);
    end
end
end

function [bp]=bpclen(A, v1, v2)
a0=besselj(0,v1); j01=a0;
b0=besselj(0,v2); j02=b0;
a1=besselj(1,v1); j11=a1;
b1=besselj(1,v2); j12=b1;
bp=A(1)*(a0.*b0);
j=2; r=1; tol=eps;
while(j<=length(A) && any(r(:)>tol))
    [a0, a1]=deal(a1, jn(j,v1,j01,j11));
    [b0, b1]=deal(b1, jn(j,v2,j02,j12));
    delta=A(j)*(a0.*b0);
    bp=bp+delta;
    r=abs(delta)./(abs(bp)+tol);
    j=j+1;
end
end

function [y]=jn(n, x, j0, j1)
y=(n>1)*ones(size(x)); 
yy=zeros(size(x));
for k=n-1:-1:1
    temp=y;
    y=2*k*y./x-yy;
    yy=temp;
end
y=y.*j1-yy.*j0;
end