function [U,H,N,Q,rr,th] = pnlse2(m,n,L,lam,pot)
% Propagator for non-linear Schrodinger equation, polar coordinates
% m Gauss-Legendre-Lobatto nodes
% n Fourier modes, must be even

mm=2*m;
[Dz,z,w]=legD(mm);
[zq,wq]=gauleg(-1,1,2*mm);
J=legC(z,zq);
D=J*Dz;

rq=abs(L*zq);
K=(1/L)*(D'*diag(wq.*rq)*D) + L*(J'*diag(wq.*rq.*(lam+pot(rq)))*J);
M=L*(J'*diag(wq./rq)*J);
B=L*(J'*diag(wq.*rq)*J);


% Symmetry BC
p1=m:-1:1;
p2=m+1:mm;
K=K(p1,p1)+K(p1,p2)+K(p2,p1)+K(p2,p2);
M=M(p1,p1)+M(p1,p2)+M(p2,p1)+M(p2,p2);
B=B(p1,p1)+B(p1,p2)+B(p2,p1)+B(p2,p2);
J=J(mm+1:2*mm,m+1:2*m);
jac=L*wq(mm+1:2*mm).*rq(mm+1:2*mm);


[rr,th]=ndgrid(L*z(m+1:mm),(2*pi/n)*(0:n-1));

% Dirichlet BC
rd=[];
kd=setdiff(1:m,rd);

% Fast diagonalization
function [V,E]=fdm(A,B,kd)
    E=zeros(m,1);
    V=eye(m);
    [V(kd,kd),E(kd)]=eig(A(kd,kd),B(kd,kd),'vector');
    s=sqrt(diag(V'*B*V));
    V=V/diag(s);
end

% Separation of variables
kt=fftshift((-n/2):(n/2-1));

lt=kt.^2;
V=zeros(m,m,n);
E=zeros(m,n);
for k=1:n
    [V(:,:,k),E(:,k)]=fdm(K+lt(k)*M,B,kd);
end


function u=prop(dt,u)
% u = V*exp(1i*dt*E)*V'*M*u
    u=fft(B*u,[],2);
    for j=1:n
       u(:,j)=V(:,:,j)'*u(:,j); 
    end
    u=exp(-1i*dt*E).*u;
    for j=1:n
       u(:,j)=V(:,:,j)*u(:,j); 
    end
    u=ifft(u,[],2);
end


function E=energy(v,u)
% E = v'*H*u
    hu=K*u+ifft(fft(M*u,[],2)*diag(lt),[],2);
    E=v(:)'*hu(:)/(2*n);
end

function N=numop(v,u)
% E = v'*M*u
    bu=B*u;
    N=v(:)'*bu(:)/n;
end

function E=expval(v,q,u)
    ju=J*u;
    jv=J*v;
    jq=J*q;
    ju=(diag(jac)*jq).*ju;
    E=jv(:)'*ju(:)/n;
end

U=@prop;
H=@energy;
N=@numop;
Q=@expval;
end

