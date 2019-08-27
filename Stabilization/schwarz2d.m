function [A,B] = schwarz2d(n,no,nu,vx,vy,bc)
CFL=0.3;

[Dhat,xhat,what]=legD(n);
Qsvv=svvker(xhat);
Bhat=diag(what);
Ahat=Dhat'*Bhat*Dhat;
Chat=Bhat*Dhat;
Asvv=Dhat'*Qsvv*Dhat;
F=bubfilt(xhat);
HPF=Bhat*(eye(n)-F);

peh=hypot(vx,vy)*max(diff(xhat))/(2*nu);
dx=min(diff(xhat));
dt=CFL*dx;
idt=1/dt;
nuh=0/(n-1);
%nuh=dx;

Astab=nu*Ahat+nuh*Asvv+idt*HPF;

j1=1:no+1;
j0=j1(end)+1:j1(end)+n;
j2=j0(end)+1:j0(end)+1+no;

nxb=n+2*(no+1);
J=zeros(n,nxb);
J(:,j0)=eye(n);
J(:,j1)=J(:,j0(n-no:n));
J(:,j2)=J(:,j0(1:no+1));

JX=J;
JY=J;

p=1.00;
if(peh>1)
if(vx>0)
    JX(no+2:end,end)=p;
elseif(vx<0)
    JX(1:n-no-1,1)=p;
elseif(vx==0&&vy~=0)
    JX(no+2:end,end)=p;
    JX(1:n-no-1,1)=p;    
end
if(vy>0)
    JY(no+2:end,end)=p;
elseif(vy<0)
    JY(1:n-no-1,1)=p;
elseif(vy==0&&vx~=0)
    JY(no+2:end,end)=p;
    JY(1:n-no-1,1)=p;
end
end



A=zeros(nxb,nxb,2);
A(:,:,1)=J'*(Astab+vx*Chat)*JX;
A(:,:,2)=J'*Bhat*JX;
A=schwarz1d(n,no,A,bc(1),bc(2));

B=zeros(nxb,nxb,2);
B(:,:,1)=J'*(Astab+vy*Chat)*JY;
B(:,:,2)=J'*Bhat*JY;
B=schwarz1d(n,no,B,bc(3),bc(4));

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

