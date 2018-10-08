function [lambda,mu,q] = qhoelliptical(a,b,m,n,k,mode,omega)
% Solves Helmoholtz in a elliptical domain as a two-parameter eigenproblem.

% Semi-focal distance
f=sqrt(a^2-b^2);
q=omega*f^2;

% Differential operators and Jacobian
% (K1 + lam*M1 + mu*B1)*u1 = 0
% (K2 + lam*M2 + mu*B2)*u2 = 0

xi0=acosh(a/f);
[Dx,x]=legD(m);
[xq,wx]=gauleg(-1,1,2*m);
J1=legC(x,xq);
Dx=2/xi0*Dx; x=xi0/2*(x+1);
xq=xi0/2*(xq(:)+1); wx=xi0/2*wx(:);
gx=wx.*exp(-(q/2)*cosh(2*xq));
W1=J1'*diag(gx)*J1;
M1=J1'*diag(q*cosh(2*xq).*gx)*J1;
K1=Dx'*W1*Dx;
B1=-q*W1;


[Dy,y]=legD(n);
[yq,wy]=gauleg(-1,1,2*n);
J2=legC(y,yq);
Dy=4/pi*Dy; y=pi/4*(y'+1);
yq=pi/4*(yq(:)+1); wy=pi/4*wy(:);
gy=wy.*exp(-(q/2)*cos(2*yq));
W2=J2'*diag(gy)*J2;
M2=J2'*diag(-q*cos(2*yq).*gy)*J2;
K2=Dy'*W2*Dy;
B2=q*W2;


% Boundary conditions
switch mode
    case 1 % MathieuJe(2*s,q,x)*MathieuC(2*s,q,y)
        rd1 = m; % F(xi0) = F'(0) = 0
        rd2 = []; % G'(pi/2) = G'(0) = 0
    case 2 % MathieuJe(2*s+1,q,x)*MathieuC(2*s+1,q,y)
        rd1 = m; % F(xi0) = F'(0) = 0
        rd2 = n; % G(pi/2) = G'(0) = 0
    case 3 % MathieuJo(2*s+2,q,x)*MathieuS(2*s+2,q,y)
        rd1 = [1,m]; % F(xi0) = F(0) = 0
        rd2 = [1,n]; % G(pi/2) = G(0) = 0
    case 4 % MathieuJo(2*s+1,q,x)*MathieuS(2*s+1,q,y)
        rd1 = [1,m]; % F(xi0) = F(0) = 0
        rd2 = 1; % G'(pi/2) = G(0) = 0
end

kd1=setdiff(1:m,rd1);
kd2=setdiff(1:n,rd2);



if mode==1  % in first mode A2 is singular, we shift mu to mu - 5
    K1=K1+5*M1;
    K2=K2+5*M2;
end

% Compute first k eigenmodes
V1=zeros(m,k);
V2=zeros(n,k);
[mu,lambda,V1(kd1,:),V2(kd2,:)]=twopareigs(K1(kd1,kd1),B1(kd1,kd1),M1(kd1,kd1),...
                                           K2(kd2,kd2),B2(kd2,kd2),M2(kd2,kd2),k);
mu=mu-5*(mode==1);

[lambda,id]=sort(lambda,'ascend');
mu=mu(id);
V1=V1(:,id);
V2=V2(:,id);

% Interpolation
% x=linspace(real(xi0),0,1024)';
% y=linspace(pi/2,0,1024);
% V1=interpcheb(V1,linspace(1,-1,1024),1);
% V2=interpcheb(V2,linspace(1,-1,1024),1);

% Fix periodicity
V2=flipud(V2);
y=y(end:-1:1);
y=[y+3*pi/2,y+pi,y+pi/2,y];
V2=[(-1)^(mode==3||mode==2)*flipud(V2);V2];
V2=[(-1)^(mode==3||mode==4)*flipud(V2);V2];

% Plot
figure(1); plot(x,V1);
figure(2); plot(y,V2);

xx=f*cosh(x)*cos(y);
yy=f*sinh(x)*sin(y);
rr=hypot(yy,xx);
RR=exp(-omega*rr.^2/2);
uuu=zeros(length(x),length(y),k);
for i=1:k
    uuu(:,:,i)=RR.*(V1(:,i)*V2(:,i)');
end

figure(3);
modegallery(xx(:,[end,1:end]),yy(:,[end,1:end]),uuu(:,[end,1:end],:));
xlim([-a,a]);
ylim([-b,b]);
pbaspect([a,b,hypot(a,b)/sqrt(2)]);
colormap(jet(256));
%camlight; 
shading interp;
view(2);

ord=2;
deg=0;
figure(4);
v1=real(InceC(ord,deg,q,1i*x(:)));
v2=InceC(ord,deg,q,y(:));
vv=RR.*(v1*v2');

figure(4);
surf(xx(:,[end,1:end]),yy(:,[end,1:end]),vv(:,[end,1:end]));
xlim([-a,a]);
ylim([-b,b]);
pbaspect([a,b,hypot(a,b)/sqrt(2)]);
colormap(jet(256));
%camlight; 
shading interp;
view(2);
end