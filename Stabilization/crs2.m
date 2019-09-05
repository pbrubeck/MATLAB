function [A,supg] = crs2(vert,hx,hy,nu,vx,vy)

% Rectangular domains
% Piecewise constant velocities
v=hypot(vx,vy);
h=v./max(abs(vx./hx),abs(vy./hy));
h(v==0)=hypot(hx(v==0),hy(v==0));
peh=v.*h/nu;
tau=(h./v).*(1-1./peh);
tau(peh<=1)=0;

nel=size(vert,ndims(vert));
N=2;
n=N*N;
[D1,x1]=legD(N);
V1=VandermondeLeg(x1);
M1=inv(V1*V1');
K1=D1'*M1*D1;
C1=M1*D1;

I=eye(n);
M=kron(M1,M1);
Kx=kron(M1,K1);
Ky=kron(K1,M1);
Cx=kron(M1,C1);
Cy=kron(C1,M1);

aa=zeros(n,n,nel);
ss=zeros(n,n,nel);
ia=zeros(n,n,nel);
ja=zeros(n,n,nel);
kk=reshape(vert,[],nel);
for e=1:nel
    K=(hy(e)/hx(e))*Kx+(hx(e)/hy(e))*Ky;
    C=(hy(e)*vx(e))*Cx+(hx(e)*vy(e))*Cy;
    S=I+(tau(e)/(hx(e)*hy(e)))*(C'/M);
    aa(:,:,e)=nu*K+S*C;
    ss(:,:,e)=S;
    [ia(:,:,e),ja(:,:,e)]=ndgrid(kk(:,e));
end
A=sparse(ia(:),ja(:),aa(:));

% SUPG solver
function [u]=sfun(A,r)
    sz=size(r);
    r=reshape(r,n,[]);
    for ie=1:size(r,2)
        r(:,ie)=ss(:,:,ie)*r(:,ie);
    end
    w=full(sparse(vert(:),ones(numel(vert),1),r(:)));
    z=A\w;
    u=reshape(z(vert),sz);
end
supg=@(A,x) sfun(A,x);
end

