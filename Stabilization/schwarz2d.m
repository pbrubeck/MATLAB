function [SA,SB]=schwarz2d(n,no,hx,hy,nu,vx,vy,dt,bc,ifdeal,ifneu)
% Assumes rectangular elements and constant velocity 
nxb=n+2*(no+1);
ns=nxb-2;
tol=1e-6;

% SEM hat
[Dhat,xhat,what]=legD(n);
if(ifdeal)
    V=VandermondeLeg(xhat);
    Bhat=inv(V*V');
else
    Bhat=diag(what);
end
Ahat=Dhat'*Bhat*Dhat;
Chat=Bhat*Dhat;

% Stablization
HPF=Bhat*(eye(n)-bubfilt(xhat))/dt;

% grid Peclet number
peh=(vx.^2+vy.^2)./max(abs(vx./hx),abs(vy./hy))*min(diff(xhat))./(nu);
nel=numel(peh);
peh(:)=10;

% Omega_bar basis
j1=1:no+1;
j0=j1(end)+1:j1(end)+n;
j2=j0(end)+1:j0(end)+1+no;

J=zeros(n,nxb);
J(:,j0)=eye(n);
J(:,j1)=J(:,j0(n-no:n));
J(:,j2)=J(:,j0(1:no+1));

nmat=3;
A=zeros(nxb,nxb,nmat);
B=zeros(nxb,nxb,nmat);
SA=zeros(ns,ns,nmat,nel);
SB=zeros(ns,ns,nmat,nel);

for e=1:nel
nux=nu/hx(e);
nuy=nu/hy(e);

% Neumann (constant extrapolation) for convection dominated
JX=J;
JY=J;
if(ifneu)
    JX(no+2:end,end)=(vx(e)>0 || abs(vx(e))<tol && peh(e)>1);
    JX(1:n-no-1,1)  =(vx(e)<0 || abs(vx(e))<tol && peh(e)>1);
    JY(no+2:end,end)=(vy(e)>0 || abs(vy(e))<tol && peh(e)>1);
    JY(1:n-no-1,1)  =(vy(e)<0 || abs(vy(e))<tol && peh(e)>1);
end

A(:,:,1)=J'*(nux*Ahat+vx(e)*Chat+hx(e)*HPF)*JX;
A(:,:,2)=J'*(hx(e)*Bhat)*JX;
A(:,:,3)=J'*(-dt*hx(e)*HPF)*JX;
SA(:,:,:,e)=schwarz1d(n,no,A,bc(1,e),bc(2,e));

B(:,:,1)=J'*(nuy*Ahat+vy(e)*Chat+hy(e)*HPF)*JY;
B(:,:,2)=J'*(hy(e)*Bhat)*JY;
B(:,:,3)=J'*(hy(e)*HPF)*JY;
SB(:,:,:,e)=schwarz1d(n,no,B,bc(3,e),bc(4,e));
end

return;


% Debugging
e=1;
P0=kron(SB(:,:,2,e),SA(:,:,1,e))+kron(SB(:,:,1,e),SA(:,:,2,e));
save('P0','P0');
return;

Mbar=J'*Bhat*J;
Mbar=schwarz1d(n,no,Mbar,2,2);
R=Mbar*ones(size(Mbar,1),size(Mbar,1))*Mbar';

AA=A(:,:,2,1)\A(:,:,1,1);
BB=B(:,:,2,1)\B(:,:,1,1);
RR=Mbar\R/Mbar';
uu=sylvester(AA,BB',RR);

xs=[xhat(n-no:n)-2;xhat(2:n-1);xhat(1:no+1)+2];
xs=[xhat(n-no-1:n)-2;xhat(2:n-1);xhat(1:no+2)+2];
vv=uu;
uu=zeros(length(xs),length(xs));
uu(2:end-1,2:end-1)=vv;
[xx,yy]=ndgrid(xs);
figure(1);
surf(xx,yy,uu);
end

