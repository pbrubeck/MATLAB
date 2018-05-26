function [] = poissonRhomboid(m)
% Single element Nearest Kronecker Product preconditioner

% Low rank = k
k=3;

% Angle of rhomboid = phi
phi=pi/4;

% Mapping
[xv,yv]=ndgrid([-1;1]);
zvert=xv(:)+1i*yv(:);
zvert=exp(1i*phi)*zvert;
zvert=2*real(zvert)+1i*imag(zvert);
rad=zeros(size(zvert)); rad(:,:)=inf;
map=curvedquad(zvert,rad);

rd=[1,m];
kd=2:m-1;

% Set up SEM
[D0,x0]=legD(m);
D0=D0(end:-1:1,end:-1:1);
x0=x0(end:-1:1);

[x1,w1]=gauleg(-1,1,ceil(3*m/2));
E10=legC(x0,x1);

I0=eye(m);
M0=E10'*diag(w1)*E10;
K0=D0'*M0*D0;
S0=M0*D0;

% Boundary conditions
a=[1,1];
b=[0,0];
C=diag(a)*I0(rd,:)+diag(b)*D0(rd,:);
G=eye(m);
G(rd,kd)=-C(:,kd);
G(rd,:)=C(:,rd)\G(rd,:);

M=G'*M0*G;
K=G'*K0*G;
S=G'*S0*G;

% Eigenvectors
V1=zeros(m,m);
V2=zeros(m,m);
[V1(kd,kd),LK]=eig(K(kd,kd),M(kd,kd),'vector');
V1(kd,kd)=V1(kd,kd)./sqrt(sum(V1(kd,kd).*(M(kd,kd)*V1(kd,kd)),1));
[LK,p]=sort(LK);
V1(kd,kd)=V1(kd,kd(p));
V1(rd,kd)=G(rd,kd)*V1(kd,kd);

% Metric coefficients
[jac,g11,g12,g22]=diffgeom(map,0,0);

c11= g22/sqrt(jac);
c12=-g12/sqrt(jac);
c22= g11/sqrt(jac);

% Restriction matrix
R=zeros(m,m-2);
R(kd,:)=eye(m-2);

% Matrix-free operators
function [vv]=assembly(uu)
    uu=reshape(uu,[m-2,m-2]);
    vv=R*uu*R';
end

function [vv]=restrict(uu)
    vv=R'*uu*R;
end

function [vv]=lapop(uu)
    vv = c11*(K0*uu*M0')+2*c12*(S0*uu*S0')+c22*(M0*uu*K0');
end

function [vv]=fullop(uu)
    vv = restrict(lapop(assembly(uu)));
end

function [L]=lapfullop()
    I=eye((m-2)*(m-2));
    L=zeros(size(I));
    for j=1:size(I,1)
        L(:,j)=reshape(fullop(I(:,j)),[],1);
    end
end

% Full operator
L=lapfullop();

% Singular Value Decomposition
Lhat=reshape(permute(reshape(L,[m-2,m-2,m-2,m-2]),[1,3,2,4]),[(m-2)*(m-2),(m-2)*(m-2)]);
[UU,SS,VV]=svds(Lhat,k);
UU=reshape(UU,m-2,m-2,k);
VV=reshape(VV,m-2,m-2,k);
SS=diag(SS);
Lnkp=zeros(size(L));
for i=1:k
    Lnkp = Lnkp + SS(i)*kron(UU(:,:,i),VV(:,:,i));  
end

figure(2);
imagesc(log(abs(L)));
colormap(gray(256)); colorbar();
title('L');

figure(3);
imagesc(log(abs(Lhat)));
colormap(gray(256)); colorbar();
title('Lhat');

figure(4);
imagesc(log(abs(Lnkp)));
colormap(gray(256)); colorbar();
title('Lnkp');

figure(5);
imagesc(log(abs(L-Lnkp)));
colormap(gray(256)); colorbar();
title('L-Lnkp');

[xx,yy]=ndgrid(x0);
zz=map(xx,yy);
figure(1);
mesh(real(zz),imag(zz),0*real(zz));
view(2);
axis equal;
end

