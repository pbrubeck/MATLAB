function [A,B]=schwarz2d(n,no,hx,hy,nu,vx,vy,CFL,bc)
nux=nu/hx;
nuy=nu/hy;

% SEM hat
[Dhat,xhat,what]=legD(n);
Bhat=diag(what);
V=VandermondeLeg(xhat);
Bhat=inv(V*V');

Ahat=Dhat'*Bhat*Dhat;
Chat=Bhat*Dhat;

% Stablization
F=bubfilt(xhat);
HPF=Bhat*(eye(n)-F);
Qsvv=svvker(xhat);
Asvv=Dhat'*Qsvv*Dhat;

peh=hypot(vx,vy)*hypot(hx,hy)*max(diff(xhat))/(nu);
dx=min(diff(xhat));
dt=CFL*dx;
idt=1/dt;
nuh=0/(n-1);
%nuh=dx;
Astab=nuh*Asvv+idt*HPF;

% Omega_bar basis
j1=1:no+1;
j0=j1(end)+1:j1(end)+n;
j2=j0(end)+1:j0(end)+1+no;

nxb=n+2*(no+1);
J=zeros(n,nxb);
J(:,j0)=eye(n);
J(:,j1)=J(:,j0(n-no:n));
J(:,j2)=J(:,j0(1:no+1));

% Neumann (constant extrapolation) for convection dominated
JX=J;
JY=J;

%peh=0;
if(peh>1)
if(vx>0)
    JX(no+2:end,end)=1;
elseif(vx<0)
    JX(1:n-no-1,1)=1;
elseif(vx==0&&vy~=0)
    %JX(no+2:end,end)=1;
    %JX(1:n-no-1,1)=1;    
end
if(vy>0)
    JY(no+2:end,end)=1;
elseif(vy<0)
    JY(1:n-no-1,1)=1;
elseif(vy==0&&vx~=0)
    %JY(no+2:end,end)=1;
    %JY(1:n-no-1,1)=1;
end
end

A=zeros(nxb,nxb,2);
A(:,:,1)=J'*(Astab+nux*Ahat+vx*Chat)*JX;
A(:,:,2)=J'*(hx*Bhat)*JX;
A=schwarz1d(n,no,A,bc(1),bc(2));

B=zeros(nxb,nxb,2);
B(:,:,1)=J'*(Astab+nuy*Ahat+vy*Chat)*JY;
B(:,:,2)=J'*(hy*Bhat)*JY;
B=schwarz1d(n,no,B,bc(3),bc(4));

return;


% Debugging
Mbar=J'*Bhat*J;
Mbar=schwarz1d(n,no,Mbar,2,2);
R=Mbar*ones(size(Mbar,1),size(Mbar,1))*Mbar';

AA=A(:,:,2)\A(:,:,1);
BB=B(:,:,2)\B(:,:,1);
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

