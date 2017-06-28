function [  ] = testmap( m,n )
if nargin<2
    n=m;
end
[Dx,x]=chebD(n);
[Dy,y]=chebD(m); y=y';
[xx,yy]=ndgrid(x,y);

[ex,ey]=ndgrid(linspace(-1,1,16));
Z=ex+1i*ey;
Z=Z./abs(Z).*max(abs(real(Z)),abs(imag(Z)));

em=size(Z,1)-1;
en=size(Z,2)-1;
zz=zeros(m, em, n, en);
jac=zeros(m, em, n, en);
for i=1:em
    for j=1:en
        [f,J]=mapquad(Z(i:i+1,j:j+1));
        zz(:,i,:,j)=f(xx,yy);
        jac(:,i,:,j)=J(xx,yy);
    end
end

zz=reshape(zz, em*m, en*n);
jac=reshape(jac, em*m, en*n);

figure(1);
surf(real(zz), imag(zz), jac);
axis square;
view(2);
end

