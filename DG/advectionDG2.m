function [] = advectionDG2(k,p)
k(1:2)=k; k1=k(1); k2=k(2);
p(1:2)=p; p1=p(1); p2=p(2);
m=k1*p1;
n=k2*p2;

% Multiplication times block-diagonals
blockx=@(A,X) reshape(A*reshape(X  ,p1,[]),size(X));
blocky=@(A,X) reshape(A*reshape(X.',p2,[]),fliplr(size(X))).';

% Diff matrix, nodes and quadrature
[Dx,xn,wx]=legD(p1); 
[Dy,yn,wy]=legD(p2);
yn=yn';

% Elements
xe=linspace(-1,1,k1+1)';
ye=linspace(-1,1,k2+1)';
[xxe,yye]=ndgrid(xe,ye);
Z=(xxe+1i*yye);

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
Q=[1,1;1,-1]/2;
zz=zeros(p1,k1,p2,k2);
LL=zeros(p1,k1,p2,k2);
JJ=zeros(p1,k1,p2,k2);
zx=zeros(p1,k1,p2,k2);
zy=zeros(p1,k1,p2,k2);
for i=1:k1
    for j=1:k2
        A=Q*Z(i+1:-1:i,j+1:-1:j)*Q;
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
mass=@(X) mod2nod(nod2mod(X)./LL);

% Stiffness matrix
M1=inv(U1*U1');
M2=inv(U2*U2');
K1=(U1*U1')*diag(wx)*Dx;
K2=(U2*U2')*diag(wy)*Dy;

stiff=@(X) blocky(M2,blockx(M1,imag(conj(zx).*blockx(K1,X)-conj(zy).*blocky(K2,X))));

% Lift operator
lift=@(X) real(blocky(M2,conj(-zx).*upwindx(X))-blockx(M1,conj(zy).*upwindy(X)));

% Advection velocity
c=1+1i;
s=1;

% Upwind flux
function [F]=upwindx(q)
    F=zeros(size(q));
    for ix=1:n
        qq=reshape(q(:,ix),p1,[]);
        q0=(qq(end,[end,1:end-1])+qq(1,:))/2;
        qj=(qq(end,[end,1:end-1])-qq(1,:))/2;
        G=zeros(size(qq));
        f=c*q0+s*abs(real(c))*qj;
        G(end,[end,1:end-1])=f;
        G(1,:)=-f;
        F(:,ix)=reshape(G,m,[]);
    end
end

function [F]=upwindy(q)
    F=zeros(size(q));
    for ix=1:m
        qq=reshape(q(ix,:).',p2,[]);
        q0=(qq(end,[end,1:end-1])+qq(1,:))/2;
        qj=(qq(end,[end,1:end-1])-qq(1,:))/2;
        G=zeros(size(qq));
        f=c*q0+s*1i*abs(imag(c))*qj;
        G(end,[end,1:end-1])=f;
        G(1,:)=-f;
        F(ix,:)=reshape(G,[],n);
    end
end

% Flux
function [F]=linearflux(u)
    F=c*u;
end

xx=real(zz);
yy=imag(zz);
u=0.1*sin(2*pi*xx).*(exp(-20*(yy-0.5).^2)+exp(-20*(yy+0.5).^2));

% To do: fix signs!
H=kron((-1).^(1:k1)'*(-1).^(1:k2),ones(p1,p2));
u=H.*u;

% Plot
figure(1);
h0=surf(real(zz), imag(zz), H.*u(:,:,1));
shading interp; %camlight;
colormap(jet(256)); colorbar;
axis square manual; caxis manual;
view(2);
drawnow;

function du=partialT(u)
    F=linearflux(u);
    du=zeros(size(u));
    for l=1:size(u,3)
        du(:,:,l)=mass(stiff(F(:,:,l))-lift(u(:,:,l)));
    end
end

function u=solveRK4(u,dt)
    rk1=dt*partialT(u);
    rk2=dt*partialT(u+rk1/2);
    rk3=dt*partialT(u+rk2/2);
    rk4=dt*partialT(u+rk3);
    u=u+(rk1+2*rk2+2*rk3+rk4)/6;
end



dt=0.5*min(1/m,1/n)/abs(c);
T=3;
nframes=ceil(T/dt);
for i=1:nframes
    u=solveRK4(u,dt);
    set(h0,'ZData',H.*u(:,:,1));
    drawnow;
end

end
