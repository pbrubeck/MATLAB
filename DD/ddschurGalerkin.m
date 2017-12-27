function [S,H,gf,K1,K2,M1,M2,x,y] = ddschurGalerkin(adjx,adjy,m,n,a,b)
% Computes the Schur complement for a separable operator kron(M2,K1)+kron(K2,M1)
rd1=[1,m];
rd2=[1,n];
kd1=2:m-1;
kd2=2:n-1;
east =rd1(1);
west =rd1(2);
north=rd2(1);
south=rd2(2);

[K1,M1,E1,V1,L1,x]=GalerkinGLL(m,a(1,:),b(1,:));
[K2,M2,E2,V2,L2,y]=GalerkinGLL(n,a(2,:),b(2,:));
[Lx,Ly]=ndgrid(L1,L2);
LL=1./(Lx+Ly);
LL((Lx==0)|(Ly==0))=0;
function uu=greenF(F,b1,b2)
    uu=zeros(m,n);
    uu(rd1,kd2)=b1;
    uu(kd1,rd2)=b2;
    rhs=M1*F*M2'-(K1*uu*M2'+M1*uu*K2');
    uu(kd1,kd2)=V1(kd1,kd1)*((V1(kd1,kd1)'*rhs(kd1,kd2)*V2(kd2,kd2)).*LL)*V2(kd2,kd2)';
    uu=E1*uu*E2';
end
gf=@(F,b1,b2) greenF(F,b1,b2);

% Schur complement matrix

% Building blocks

Q1=M1(kd1,kd1)*V1(kd1,kd1);
Q2=M2(kd2,kd2)*V2(kd2,kd2);

% S11
p1=[east,east,west,west];
p2=[east,west,east,west];
VTML=V1(kd1,kd1)'*M1(p1,kd1)';
VTKL=V1(kd1,kd1)'*K1(p1,kd1)';
VTMR=V1(kd1,kd1)'*M1(kd1,p2) ;
VTKR=V1(kd1,kd1)'*K1(kd1,p2) ;
MM=VTML.*VTMR;
MK=VTML.*VTKR;
KM=VTKL.*VTMR;
KK=VTKL.*VTKR;
DXX=diag(L2.^2)*(LL'*MM)+diag(L2)*(LL'*(MK+KM))+(LL'*KK);
EE=Q2*diag(DXX(:,1))*Q2';
EW=Q2*diag(DXX(:,2))*Q2';
WE=Q2*diag(DXX(:,3))*Q2';
WW=Q2*diag(DXX(:,4))*Q2';

% S22
p1=[north,north,south,south];
p2=[north,south,north,south];
VTML=V2(kd2,kd2)'*M2(p1,kd2)';
VTKL=V2(kd2,kd2)'*K2(p1,kd2)';
VTMR=V2(kd2,kd2)'*M2(kd2,p2) ;
VTKR=V2(kd2,kd2)'*K2(kd2,p2) ;
MM=VTML.*VTMR;
MK=VTML.*VTKR;
KM=VTKL.*VTMR;
KK=VTKL.*VTKR;
DYY=diag(L1.^2)*(LL*MM)+diag(L1)*(LL*(MK+KM))+(LL*KK);
NN=Q1*diag(DYY(:,1))*Q1';
NS=Q1*diag(DYY(:,2))*Q1';
SN=Q1*diag(DYY(:,3))*Q1';
SS=Q1*diag(DYY(:,4))*Q1';

% S12
p1=[east, east, west, west ];
p2=[north,south,north,south];
VTML=V1(kd1,kd1)'*M1(p1,kd1)';
VTKL=V1(kd1,kd1)'*K1(p1,kd1)';
VTMR=V2(kd2,kd2)'*M2(kd2,p2) ;
VTKR=V2(kd2,kd2)'*K2(kd2,p2) ;
U=[L2.*VTMR, VTMR, L2.*VTKR, VTKR];
V=[L1.*VTML, L1.*VTKL, VTML, VTKL];
EN=Q1*((U(:,1:4:end)*V(:,1:4:end)').*LL')*Q2';
ES=Q1*((U(:,2:4:end)*V(:,2:4:end)').*LL')*Q2';
WN=Q1*((U(:,3:4:end)*V(:,3:4:end)').*LL')*Q2';
WS=Q1*((U(:,4:4:end)*V(:,4:4:end)').*LL')*Q2';

% S21
p1=[north,north,south,south];
p2=[east, west, east, west ];
VTML=V2(kd2,kd2)'*M2(p1,kd2)';
VTKL=V2(kd2,kd2)'*K2(p1,kd2)';
VTMR=V1(kd1,kd1)'*M1(kd1,p2) ;
VTKR=V1(kd1,kd1)'*K1(kd1,p2) ;
U=[L1.*VTMR, VTMR, L1.*VTKR, VTKR];
V=[L2.*VTML, L2.*VTKL, VTML, VTKL];
NE=Q2*((U(:,1:4:end)*V(:,1:4:end)').*LL)*Q1';
NW=Q2*((U(:,2:4:end)*V(:,2:4:end)').*LL)*Q1';
SE=Q2*((U(:,3:4:end)*V(:,3:4:end)').*LL)*Q1';
SW=Q2*((U(:,4:4:end)*V(:,4:4:end)').*LL)*Q1';

% Assembly

% S11
[x1,y1]=ndgrid(adjx(:,1), adjx(:,1));
[x2,y2]=ndgrid(adjx(:,2), adjx(:,2));

S11=sparse((n-2)*size(adjx,1),(n-2)*size(adjx,1));
S11=S11+kron(sparse(x2==y2), -EE+M1(east,east)*K2(kd2,kd2)+K1(east,east)*M2(kd2,kd2));
S11=S11+kron(sparse(x2==y1), -EW+M1(east,west)*K2(kd2,kd2)+K1(east,west)*M2(kd2,kd2));
S11=S11+kron(sparse(x1==y2), -WE+M1(west,east)*K2(kd2,kd2)+K1(west,east)*M2(kd2,kd2));
S11=S11+kron(sparse(x1==y1), -WW+M1(west,west)*K2(kd2,kd2)+K1(west,west)*M2(kd2,kd2));
            
% S22
[x1,y1]=ndgrid(adjy(:,1), adjy(:,1));
[x2,y2]=ndgrid(adjy(:,2), adjy(:,2));

S22=sparse((m-2)*size(adjy,1),(m-2)*size(adjy,1));
S22=S22+kron(sparse(x2==y2), -NN+M2(north,north)*K1(kd1,kd1)+K2(north,north)*M1(kd1,kd1));
S22=S22+kron(sparse(x2==y1), -NS+M2(north,south)*K1(kd1,kd1)+K2(north,south)*M1(kd1,kd1));
S22=S22+kron(sparse(x1==y2), -SN+M2(south,north)*K1(kd1,kd1)+K2(south,north)*M1(kd1,kd1));
S22=S22+kron(sparse(x1==y1), -SS+M2(south,south)*K1(kd1,kd1)+K2(south,south)*M1(kd1,kd1));

            
% S12
[x1,y1]=ndgrid(adjx(:,1), adjy(:,1));
[x2,y2]=ndgrid(adjx(:,2), adjy(:,2));

S12=sparse((n-2)*size(adjx,1),(m-2)*size(adjy,1));
S12=S12+kron(sparse(x2==y2), -EN+M2(kd2,north)*K1(east,kd1)+K2(kd2,north)*M1(east,kd1));
S12=S12+kron(sparse(x2==y1), -ES+M2(kd2,south)*K1(east,kd1)+K2(kd2,south)*M1(east,kd1));
S12=S12+kron(sparse(x1==y2), -WN+M2(kd2,north)*K1(west,kd1)+K2(kd2,north)*M1(west,kd1));
S12=S12+kron(sparse(x1==y1), -WS+M2(kd2,south)*K1(west,kd1)+K2(kd2,south)*M1(west,kd1));

% S21
[y1,x1]=ndgrid(adjy(:,1), adjx(:,1));
[y2,x2]=ndgrid(adjy(:,2), adjx(:,2));

S21=sparse((m-2)*size(adjy,1),(n-2)*size(adjx,1));
S21=S21+kron(sparse(y2==x2), -NE+M1(kd1,east)*K2(north,kd2)+K1(kd1,east)*M2(north,kd2));
S21=S21+kron(sparse(y2==x1), -NW+M1(kd1,west)*K2(north,kd2)+K1(kd1,west)*M2(north,kd2));
S21=S21+kron(sparse(y1==x2), -SE+M1(kd1,east)*K2(south,kd2)+K1(kd1,east)*M2(south,kd2));
S21=S21+kron(sparse(y1==x1), -SW+M1(kd1,west)*K2(south,kd2)+K1(kd1,west)*M2(south,kd2));

id1=1:(n-2)*size(adjx,1);
id2=(n-2)*size(adjx,1)+(1:(m-2)*size(adjy,1));

sz=(n-2)*size(adjx,1)+(m-2)*size(adjy,1);
S=sparse(sz,sz);
S(id1,id1)=S11;
S(id1,id2)=S12;
S(id2,id1)=S21;
S(id2,id2)=S22;

H=0;
end


function [K, M, E, V, L, x]=GalerkinGLL(m,a,b)
% Stiffness and mass matrices for a Gauss-Legendre-Lobatto nodal basis over
% the interval [-1,1] with Robin BCs specified in a and b.
I=eye(m);                           % Identity
kd=2:m-1;                           % kept DOFs
rd=[1,m];                           % removed DOFs
[D,x,w]=legD(m);                    
D=D(end:-1:1,end:-1:1);             % Differentiation matrix
x=x(end:-1:1);                      % Collocation nodes
w=w(end:-1:1);                      % Quadrature
C=diag(a)*I(rd,:)+diag(b)*D(rd,:);	% Constraint operator

% Constrained basis
E=eye(m);
E(rd,kd)=-C(:,kd);
E(rd,:)=C(:,rd)\E(rd,:);

% Primed basis
DE=D*E;

% Mass matrix
V=VandermondeLeg(x);
Minv=(V*V');
M=E'*(Minv\E);

% Stiffness matrix
K=DE'*diag(w)*DE-E(rd,:)'*diag([-1,1])*DE(rd,:);

% Eigenfunctions
V=zeros(m);
[V(kd,kd),L]=eig(K(kd,kd), M(kd,kd), 'vector');
V(kd,kd)=normc(V(kd,kd), M(kd,kd));
V(rd,kd)=E(rd,kd)*V(kd,kd);
end