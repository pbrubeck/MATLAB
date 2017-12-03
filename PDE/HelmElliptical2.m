function [] = HelmElliptical2(a, b, m, n, k, mode)
% Solves Helmoholtz in a elliptical domain as a two-parameter eigenproblem.

% Semi-focal distance
f=sqrt(a^2-b^2);

% Kept and removed degrees of freedom
kd1=2:m;
kd2=1:n;
rd1=1;

% Differential operators and Jacobian
[Dx,x]=chebD(2*m); 
xi0=acosh(a/f); x=xi0*x;
Dx=1/xi0*Dx; Dxx=Dx*Dx;
x=x(1:m); s=1;
Dxx=Dxx(1:m,1:m)+s*Dxx(1:m,end:-1:m+1);
Dx=Dx(1:m,1:m)+s*Dx(1:m,end:-1:m+1);

Jx=diag(f^2/2*cosh(2*x));
E1=eye(m);

[Dyy,y]=fourD2(n);
Jy=-diag(f^2/2*cos(2*y));
E2=eye(n);

J=bsxfun(@plus, f^2/2*cosh(2*x), -f^2/2*cos(2*y));
lap=@(uu) (Dxx*uu+uu*Dyy');

% Boundary conditions
bc1=[1,0];
BC1=diag(bc1(1))*E1(rd1,:)+diag(bc1(2))*Dx(rd1,:);

% Give-back matrix
G1=-BC1(:,rd1)\BC1(:,kd1);

% Schur complements
A1=Dxx(kd1,kd1)+Dxx(kd1,rd1)*G1;
B1=Jx(kd1,kd1)+Jx(kd1,rd1)*G1;
C1=E1(kd1,kd1)+E1(kd1,rd1)*G1;

A2=Dyy;
B2=Jy;
C2=-E2;

% in first mode A2 is singular, we shift mu to mu - 5
A1=A1+5*C1;
A2=A2+5*C2;


% Compute first k eigenmodes
V1=zeros(m,k);
[mu,lambda,V1(kd1,:),V2]=twopareigs(A1,C1,B1,A2,C2,B2,k);
mu=mu-5;

[lambda,id]=sort(lambda,'descend');
mu=mu(id);
V1=V1(:,id);
V2=V2(:,id);

% Retrieve boundary values
V1(rd1,:)=G1*V1(kd1,:);

q=-f^2/4*lambda;
display(lambda);
display(mu);
display(q);

% Interpolation
% x=linspace(real(xi0),0,1024)';
% V1=interpcheb(V1,linspace(1,-1,1024),1);


% Plot
figure(1); plot(x,V1);
figure(2); plot(y,V2(:,end));

xx=f*cosh(x)*cos(y);
yy=f*sinh(x)*sin(y);
uuu=zeros(m,n,k);
for i=1:k
    uuu(:,:,i)=V1(:,i)*V2(:,i)';
end

figure(3);
modegallery(xx(:,[end,1:end]),yy(:,[end,1:end]),uuu(:,[end,1:end],:));
colormap(jet(256));
camlight; shading interp;
view(2);
end