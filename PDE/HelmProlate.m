function [] = HelmProlate(f, a, m, n, k, mode)
% Solves Helmoholtz in a elliptical domain as a two-parameter eigenproblem.

% Semi-focal distance
a0=acosh(a/f);

% Kept and removed degrees of freedom
kd1=2:m-1;
kd2=2:n-1;
rd1=[1,m];
rd2=[1,n];

% Differential operators and Jacobian
[Dx,x]=chebD(m);
x=a0/2*(x+1);
Dx=2/a0*Dx; Dxx=(Dx+diag(coth(x)))*Dx;
Jx=diag(a0^2*sinh(x).^2);
E1=eye(m);

[Dy,y]=chebD(n);
y=pi/4*(y'+1);
Dy=4/pi*Dy; Dyy=(Dy+diag(cot(y)))*Dy;
Jy=diag(a0^2*sin(y).^2);
E2=eye(n);

J=bsxfun(@plus, a0^2/2*cosh(2*x), -a0^2/2*cos(2*y));
lap=@(uu) (Dxx*uu+uu*Dyy');

% Boundary conditions
switch mode
    case 1 % Even
        bc1 = [1 0; 0 1]; % F(a0) = F'(0) = 0
        bc2 = [0 1; 0 1]; % G'(pi/2) = G'(0) = 0
    case 2 % Odd
        bc1 = [1 0; 0 1]; % F(a0) = F'(0) = 0
        bc2 = [1 0; 0 1]; % G(pi/2) = G'(0) = 0
end
BC1=diag(bc1(:,1))*E1(rd1,:)+diag(bc1(:,2))*Dx(rd1,:);
BC2=diag(bc2(:,1))*E2(rd2,:)+diag(bc2(:,2))*Dy(rd2,:);

% Give-back matrix
G1=-BC1(:,rd1)\BC1(:,kd1);
G2=-BC2(:,rd2)\BC2(:,kd2);

% Schur complements
A1=Dxx(kd1,kd1)+Dxx(kd1,rd1)*G1;
B1=Jx(kd1,kd1)+Jx(kd1,rd1)*G1;
C1=E1(kd1,kd1)+E1(kd1,rd1)*G1;

A2=Dyy(kd2,kd2)+Dyy(kd2,rd2)*G2;
B2=Jy(kd2,kd2)+Jy(kd2,rd2)*G2;
C2=-E2(kd2,kd2)-E2(kd2,rd2)*G2;

if mode==1  % in first mode A2 is singular, we shift mu to mu - 5
    A1=A1+5*C1;
    A2=A2+5*C2;
end

% Compute first k eigenmodes
V1=zeros(m,k);
V2=zeros(n,k);
[mu,lambda,V1(kd1,:),V2(kd2,:)]=twopareigs(A1,C1,B1,A2,C2,B2,k);
mu=mu-5*(mode==1);

[lambda,id]=sort(lambda,'descend');
mu=mu(id);
V1=V1(:,id);
V2=V2(:,id);

% Retrieve boundary values
V1(rd1,:)=G1*V1(kd1,:);
V2(rd2,:)=G2*V2(kd2,:);

display(lambda);
display(mu);

% Interpolation
x=linspace(a0,0,1024)';
y=linspace(pi/2,0,1024);
V1=interpcheb(V1,linspace(1,-1,1024),1);
V2=interpcheb(V2,linspace(1,-1,1024),1);

% Fix periodicity
y=[y+3*pi/2,y+pi,y+pi/2,y];
V2=[(-1)^(mode==2)*flipud(V2);V2];
V2=[               flipud(V2);V2];

% Plot
figure(1); plot(x,V1);
figure(2); plot(y,V2);

zz=f*cosh(x)*cos(y);
xx=f*sinh(x)*sin(y);
uu=V1(:,end)*V2(:,end)';

figure(3);
surf(xx,zz,uu);
colormap(jet(256));
camlight; shading interp;
view(2); axis square;
end