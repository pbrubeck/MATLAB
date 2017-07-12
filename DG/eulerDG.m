function [] = eulerDG(k,p)
k(1:2)=k; k1=k(1); k2=k(2);
p(1:2)=p; p1=p(1); p2=p(2);
m=k1*p1;
n=k2*p2;

% Multiplication times block-diagonals
blockx=@(A,X) reshape(A*reshape(X  ,p1,[]),size(X));
blocky=@(A,X) reshape(A*reshape(X.',p2,[]),fliplr(size(X))).';

% Diff matrix and Lobatto quadrature
[Dx,xn,wx]=legD(p1);
[Dy,yn,wy]=legD(p2);
yn=yn';

% Elements
xe=linspace(-1,1,k1+1)';
ye=linspace(-1,1,k2+1)';
[xxe,yye]=ndgrid(xe,ye);
Z=xxe+(1i)*yye;

% Legendre nodes
d=@(j) j./sqrt(4*j.^2-1);
[L1,W1]=trideigs(zeros(p1,1),d(1:p2-1));
[L2,W2]=trideigs(zeros(p2,1),d(1:p2-1));

% Nodal to modal transform
U1=VandermondeLeg(xn)*W1;
U2=VandermondeLeg(yn)*W2;
nod2mod=@(X) blocky(U2',blockx(U1',X));
mod2nod=@(X) blocky(U2 ,blockx(U1 ,X));

% Jacobian
Q=[1,1;-1,1]/2;
zz=zeros(p1,k1,p2,k2);
LL=zeros(p1,k1,p2,k2);
JJ=zeros(p1,k1,p2,k2);
zx=zeros(p1,k1,p2,k2);
zy=zeros(p1,k1,p2,k2);
for i=1:k1
    for j=1:k2
        A=Q*Z(i:i+1,j:j+1)*Q';
        J=imag(conj(A(:,2))*A(2,:));
        zz(:,i,:,j)=A(1,1)+bsxfun(@plus, A(2,1)*xn, A(1,2)*yn)+A(2,2)*(xn*yn);
        LL(:,i,:,j)=J(1,1)+bsxfun(@plus, J(2,1)*L1, J(1,2)*L2');
        JJ(:,i,:,j)=J(1,1)+bsxfun(@plus, J(2,1)*xn, J(1,2)*yn );
        zx(:,i,:,j)=repmat(A(1,2)+A(2,2)*xn,[1,p2]);
        zy(:,i,:,j)=repmat(A(2,1)+A(2,2)*yn,[p1,1]);
    end
end
zz=reshape(zz,m,n);
LL=reshape(LL,m,n);
JJ=reshape(JJ,m,n);
zx=reshape(zx,m,n);
zy=reshape(zy,m,n);


% Differential operators
grad=@(X) 1i*(zx.*blockx(Dx,X)-zy.*blocky(Dy,X))./JJ;

% Inverse mass matrix
imass=@(X) mod2nod(nod2mod(X)./LL);

% Stiffness and lift matrices
M1=inv(U1*U1');
M2=inv(U2*U2');
K1=Dx'*diag(wx)*(U1*U1');
K2=Dy'*diag(wy)*(U2*U2');

stiff=@(X) imag(blockx(K1,conj(zx).*X)-blocky(K2,conj(zy).*X));
stiff=@(X) stiff(blockx(M1,blocky(M2,X)));
lift =@(X) imag(blocky(M2,conj(zx).*upwindx(X))-blockx(M1,conj(zy).*upwindy(X)));


% Flux
function [F]=flux(q)
    rho=q(:,:,1);
    rhou=q(:,:,2);
    rhov=q(:,:,3);
    E=q(:,:,4);
    P=(2/5)*(E-(rhou.^2+rhov.^2)./(2*rho));
    w=(rhou+1i*rhov)./rho;
    F(:,:,1)=rho.*w;
    F(:,:,2)=rhou.*w+P;
    F(:,:,3)=rhov.*w+1i*P;
    F(:,:,4)=(E+P).*w;
end

% Upwind flux
function [F]=upwindx(u)
    umid=(u([end,p1:p1:end-1],:,:)+u(1:p1:end,:,:))/2;
    ujmp=(u([end,p1:p1:end-1],:,:)-u(1:p1:end,:,:))/2;
    F=zeros(size(u));
    F([end,p1:p1:end-1],:,:)=flux(umid);
    F(1:p1:end,:,:)=-F([end,p1:p1:end-1],:,:);
end

function [F]=upwindy(u)
    umid=(u(:,[end,p2:p2:end-1],:)+u(:,1:p2:end,:))/2;
    ujmp=(u(:,[end,p2:p2:end-1],:)-u(:,1:p2:end,:))/2;
    F=zeros(size(u));
    F(:,[end,p2:p2:end-1],:)=flux(umid);
    F(:,1:p2:end,:)=-F(:,[end,p2:p2:end-1],:);
end

function [F]=numflux(X)
    Fx=upwindx(X);
    Fy=upwindy(X);
    F=zeros(size(X));
    for l=1:size(X,3)
        F(:,:,l)=imag(blocky(M2,conj(zx).*Fx(:,:,l))-blockx(M1,conj(zy).*Fy(:,:,l)));
    end
end


% Initial condition
xi=@(x) 1/4*(erf(100*(0.5+x))+1).*(erf(100*(0.5-x))+1);

xx=real(zz);
yy=imag(zz);
rho0=xi(yy)+1;
u0=xi(yy)-0.5;
v0=0.1*sin(2*pi*xx).*(exp(-20*(yy-0.5).^2)+exp(-20*(yy+0.5).^2));
P0=2.5;

u(:,:,1)=rho0;
u(:,:,2)=rho0.*u0;
u(:,:,3)=rho0.*v0;
u(:,:,4)=5/2*P0+rho0.*(u0.^2+v0.^2)/2;

% Time evolution
function du=partialT(u)
    F=flux(u);
    nF=numflux(u);
    du=zeros(size(u));
    for l=1:size(u,3)
        du(:,:,l)=imass(stiff(F(:,:,l))-nF(:,:,l));
    end
end

function u=solveRK4(u,dt)
    rk1=dt*partialT(u);
    rk2=dt*partialT(u+rk1/2);
    rk3=dt*partialT(u+rk2/2);
    rk4=dt*partialT(u+rk3);
    u=u+(rk1+2*rk2+2*rk3+rk4)/6;
end

% Plot
figure(1);
h0=surf(real(zz), imag(zz), u(:,:,1));
colormap(jet(256)); colorbar;
axis square;
shading interp; %camlight;
view(2);
drawnow;

h=min(min(sqrt(abs(JJ))))*min(1-cos(pi./(p-1)));
dt=0.1*h;
T=3; time=0;
nframes=ceil(T/dt);
for i=1:nframes
    time=time+dt;
    u=solveRK4(u,dt);
    set(h0,'ZData',u(:,:,1));
    title(sprintf('t = %f',time));
    drawnow;
end
end