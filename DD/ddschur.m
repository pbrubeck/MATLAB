function [S,H,V1,V2,LL] = ddschur(adjx,adjy,A1,A2,Dx,Dy,C1,C2)
% Computes the Schur complement for a separable operator kron(I,A1)+kron(A2,I)
% Enforces continuity in the first derivative.
m=size(Dx,1);
n=size(Dy,1);
rd1=[1,m];
rd2=[1,n];
kd1=2:m-1;
kd2=2:n-1;
east =rd1(1);
west =rd1(2);
north=rd2(1);
south=rd2(2);

% Eigenfunctions
function [L,S,G]=eigenfunctions(A,C,m,rd,kd)
    % Give-back matrix
    G=-C(:,rd)\C(:,kd);
    % Eigenfunctions
    S=zeros(m,length(kd));
    [S(kd,:),L]=eig(A(kd,kd)+A(kd,rd)*G,'vector');
    S(rd,:)=G*S(kd,:);
end
[L1,V1]=eigenfunctions(A1,C1,m,rd1,kd1);
[L2,V2]=eigenfunctions(A2,C2,n,rd2,kd2);
[Lx,Ly]=ndgrid(L1,L2); LL=1./(Lx+Ly);

% Schur complement matrix

% S11
[x1,y1]=ndgrid(adjx(:,1), adjx(:,1));
[x2,y2]=ndgrid(adjx(:,2), adjx(:,2));

EE = Dx(east,kd1)*V1(kd1,:)*diag(V1(kd1,:)\A1(kd1,east))*LL;
EW = Dx(east,kd1)*V1(kd1,:)*diag(V1(kd1,:)\A1(kd1,west))*LL;
WE =-Dx(west,kd1)*V1(kd1,:)*diag(V1(kd1,:)\A1(kd1,east))*LL;
WW =-Dx(west,kd1)*V1(kd1,:)*diag(V1(kd1,:)\A1(kd1,west))*LL;

S11=sparse((m-2)*size(adjx,1),(m-2)*size(adjx,1));
S11=S11+kron(sparse(x2==y2), V1(kd1,:)*diag(-Dx(east,east)+EE)/V1(kd1,:));
S11=S11+kron(sparse(x2==y1), V1(kd1,:)*diag(-Dx(east,west)+EW)/V1(kd1,:));
S11=S11+kron(sparse(x1==y2), V1(kd1,:)*diag( Dx(west,east)+WE)/V1(kd1,:));
S11=S11+kron(sparse(x1==y1), V1(kd1,:)*diag( Dx(west,west)+WW)/V1(kd1,:));
H11=kron(Dx(west,west)*(x1==y1) + Dx(west,east)*(x1==y2) + ...
        -Dx(east,west)*(x2==y1) - Dx(east,east)*(x2==y2), speye(m-2));

% S22
[x1,y1]=ndgrid(adjy(:,1), adjy(:,1));
[x2,y2]=ndgrid(adjy(:,2), adjy(:,2));

NN = Dy(north,kd2)*V2(kd2,:)*diag(V2(kd2,:)\A2(kd2,north))*LL.';
NS = Dy(north,kd2)*V2(kd2,:)*diag(V2(kd2,:)\A2(kd2,south))*LL.';
SN =-Dy(south,kd2)*V2(kd2,:)*diag(V2(kd2,:)\A2(kd2,north))*LL.';
SS =-Dy(south,kd2)*V2(kd2,:)*diag(V2(kd2,:)\A2(kd2,south))*LL.';

S22=sparse((n-2)*size(adjy,1),(n-2)*size(adjy,1));
S22=S22+kron(sparse(x2==y2), V2(kd2,:)*diag(-Dy(north,north)+NN)/V2(kd2,:));
S22=S22+kron(sparse(x2==y1), V2(kd2,:)*diag(-Dy(north,south)+NS)/V2(kd2,:));
S22=S22+kron(sparse(x1==y2), V2(kd2,:)*diag( Dy(south,north)+SN)/V2(kd2,:));
S22=S22+kron(sparse(x1==y1), V2(kd2,:)*diag( Dy(south,south)+SS)/V2(kd2,:));
H22=kron(Dy(south,south)*(x1==y1) + Dy(south,north)*(x1==y2) + ...
        -Dy(north,south)*(x2==y1) - Dy(north,north)*(x2==y2), speye(n-2));

% S12
[x1,y1]=ndgrid(adjx(:,1), adjy(:,1));
[x2,y2]=ndgrid(adjx(:,2), adjy(:,2));

EN = ((V2(kd2,:)\A2(kd2,north))*( Dx(east,kd1)*V1(kd1,:))).*LL;
ES = ((V2(kd2,:)\A2(kd2,south))*( Dx(east,kd1)*V1(kd1,:))).*LL;
WN = ((V2(kd2,:)\A2(kd2,north))*(-Dx(west,kd1)*V1(kd1,:))).*LL;
WS = ((V2(kd2,:)\A2(kd2,south))*(-Dx(west,kd1)*V1(kd1,:))).*LL;

S12=sparse((m-2)*size(adjx,1),(n-2)*size(adjy,1));
S12=S12+kron(sparse(x2==y2), V1(kd1,:)*(EN)/V2(kd2,:));
S12=S12+kron(sparse(x2==y1), V1(kd1,:)*(ES)/V2(kd2,:));
S12=S12+kron(sparse(x1==y2), V1(kd1,:)*(WN)/V2(kd2,:));
S12=S12+kron(sparse(x1==y1), V1(kd1,:)*(WS)/V2(kd2,:));

% S21
[y1,x1]=ndgrid(adjy(:,1), adjx(:,1));
[y2,x2]=ndgrid(adjy(:,2), adjx(:,2));

NE =((V1(kd1,:)\A1(kd1,east))*( Dy(north,kd2)*V2(kd2,:))).*LL.';
NW =((V1(kd1,:)\A1(kd1,west))*( Dy(north,kd2)*V2(kd2,:))).*LL.';
SE =((V1(kd1,:)\A1(kd1,east))*(-Dy(south,kd2)*V2(kd2,:))).*LL.';
SW =((V1(kd1,:)\A1(kd1,west))*(-Dy(south,kd2)*V2(kd2,:))).*LL.';

S21=sparse((n-2)*size(adjy,1),(m-2)*size(adjx,1));
S21=S21+kron(sparse(y2==x2), V2(kd2,:)*(NE)/V1(kd1,:));
S21=S21+kron(sparse(y2==x1), V2(kd2,:)*(NW)/V1(kd1,:));
S21=S21+kron(sparse(y1==x2), V2(kd2,:)*(SE)/V1(kd1,:));
S21=S21+kron(sparse(y1==x1), V2(kd2,:)*(SW)/V1(kd1,:));

I1=1:(m-2)*size(adjx,1);
I2=(m-2)*size(adjx,1)+(1:(n-2)*size(adjy,1));

sz=(m-2)*size(adjx,1)+(n-2)*size(adjy,1);
S=sparse(sz,sz);
S(I1,I1)=S11;
S(I1,I2)=S12;
S(I2,I1)=S21;
S(I2,I2)=S22;

H=sparse(sz,sz);
H(I1,I1)=H11;
H(I2,I2)=H22;
end