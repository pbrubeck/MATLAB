function [] = eulerDG(k,p)
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
xe=linspace(1,2,k1+1)';
ye=linspace(1,2,k2+1)';
[xx,yy]=ndgrid(xe,ye);
Z=(xx+1i*yy).^2/2;

% Legendre nodes
d=@(x) x./sqrt(4*x.^2-1);
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

% Inverse mass matrix
mass=@(X) mod2nod(nod2mod(X)./LL);

% Stiffness matrix
stiffx=@(X) mass(blockx(Dx,X));
stiffy=@(X) mass(blocky(Dy,X));
grad=@(X) 1i*(zx.*blockx(Dx,X)-zy.*blocky(Dy,X))./JJ;


figure(1);
surf(real(zz), imag(zz), real(grad(imag(sin(zz)))));
view(2); axis square;
shading interp; colormap(jet(256)); camlight;


% Lax-Friedrichs numerical flux
function [F]=LaxFriedrichs(u)
    
end

function du=partialT(u)
    [F,G]=eulerflux(u);
    f=LaxFriedrichs(u);
    du=mass(stiffx(F)+stiffy(G)-f);
end




end


function u=solveRK4(u,dt)
    k1=dt*partialT(u);
    k2=dt*partialT(u+k1/2);
    k3=dt*partialT(u+k2/2);
    k4=dt*partialT(u+k3);
    u=u+(k1+2*k2+2*k3+k4)/6;
end

function [F,G]=eulerflux(q)
    rho=q(:,:,1);
    rhou=q(:,:,2);
    rhov=q(:,:,3);
    E=q(:,:,4);
    u=rhou./rho;
    v=rhov./rho;
    
    p=(2/5)*(E-(rhou.^2+rhov.^2)./(2*rho));

    F(:,:,1)=rhou;
    F(:,:,2)=rhou.*u+p;
    F(:,:,3)=rhou.*v;
    F(:,:,4)=u.*(E+p);
    
    G(:,:,1)=rhov;
    G(:,:,2)=rhov.*u;
    G(:,:,3)=rhov.*v+p;
    G(:,:,4)=v.*(E+p);
end
