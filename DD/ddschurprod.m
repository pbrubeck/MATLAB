function [S,H,V1,V2,LL] = ddschurprod(adjx,adjy,A1,A2,Dx,Dy,C1,C2)
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
function [L,V,G]=eigenfunctions(A,C,m,rd,kd)
    % Give-back matrix
    G=-C(:,rd)\C(:,kd);
    % Eigenfunctions
    V=zeros(m);
    L=zeros(m,1);
    %V(:,rd)=null(A);
    [V(kd,kd),L(kd)]=eig(A(kd,kd)+A(kd,rd)*G,'vector');
    V(rd,kd)=G*V(kd,kd);
end
[L1,V1]=eigenfunctions(A1,C1,m,rd1,kd1);
[L2,V2]=eigenfunctions(A2,C2,n,rd2,kd2);
[Lx,Ly]=ndgrid(L1,L2); 
LL=1./(Lx.*Ly);
LL((Lx==0)|(Ly==0))=0;

% Schur complement matrix

% S11
[x1,y1]=ndgrid(adjx(:,1), adjx(:,1));
[x2,y2]=ndgrid(adjx(:,2), adjx(:,2));

S11=sparse((n-2)*size(adjx,1),(n-2)*size(adjx,1));
S11=S11+kron(sparse(x2==y2), speye(n-2)*(-Dx(east,east)+Dx(east,kd1)*(A1(kd1,kd1)\A1(kd1,east))));
S11=S11+kron(sparse(x2==y1), speye(n-2)*(-Dx(east,west)+Dx(east,kd1)*(A1(kd1,kd1)\A1(kd1,west))));
S11=S11+kron(sparse(x1==y2), speye(n-2)*( Dx(west,east)-Dx(west,kd1)*(A1(kd1,kd1)\A1(kd1,east))));
S11=S11+kron(sparse(x1==y1), speye(n-2)*( Dx(west,west)-Dx(west,kd1)*(A1(kd1,kd1)\A1(kd1,west))));
H11=kron(Dx(west,west)*(x1==y1) + Dx(west,east)*(x1==y2) + ...
        -Dx(east,west)*(x2==y1) - Dx(east,east)*(x2==y2), speye(n-2)); 

% S22
[x1,y1]=ndgrid(adjy(:,1), adjy(:,1));
[x2,y2]=ndgrid(adjy(:,2), adjy(:,2));

S22=sparse((m-2)*size(adjy,1),(m-2)*size(adjy,1));
S22=S22+kron(sparse(x2==y2), speye(m-2)*(-Dy(north,north)+Dy(north,kd2)*(A2(kd2,kd2)\A2(kd2,north))));
S22=S22+kron(sparse(x2==y1), speye(m-2)*(-Dy(north,south)+Dy(north,kd2)*(A2(kd2,kd2)\A2(kd2,south))));
S22=S22+kron(sparse(x1==y2), speye(m-2)*( Dy(south,north)-Dy(south,kd2)*(A2(kd2,kd2)\A2(kd2,north))));
S22=S22+kron(sparse(x1==y1), speye(m-2)*( Dy(south,south)-Dy(south,kd2)*(A2(kd2,kd2)\A2(kd2,south))));
H22=kron(Dy(south,south)*(x1==y1) + Dy(south,north)*(x1==y2) + ...
        -Dy(north,south)*(x2==y1) - Dy(north,north)*(x2==y2), speye(m-2)); 

% S12
[x1,y1]=ndgrid(adjx(:,1), adjy(:,1));
[x2,y2]=ndgrid(adjx(:,2), adjy(:,2));

S12=sparse((n-2)*size(adjx,1),(m-2)*size(adjy,1));
S12=S12+kron(sparse(x2==y2), (A2(kd2,kd2)\A2(kd2,north))*Dx(east,kd1));
S12=S12+kron(sparse(x2==y1), (A2(kd2,kd2)\A2(kd2,south))*Dx(east,kd1));
S12=S12+kron(sparse(x1==y2),-(A2(kd2,kd2)\A2(kd2,north))*Dx(west,kd1));
S12=S12+kron(sparse(x1==y1),-(A2(kd2,kd2)\A2(kd2,south))*Dx(west,kd1));

% S21
[y1,x1]=ndgrid(adjy(:,1), adjx(:,1));
[y2,x2]=ndgrid(adjy(:,2), adjx(:,2));

S21=sparse((m-2)*size(adjy,1),(n-2)*size(adjx,1));
S21=S21+kron(sparse(y2==x2), (A1(kd1,kd1)\A1(kd1,east))*Dy(north,kd2));
S21=S21+kron(sparse(y2==x1), (A1(kd1,kd1)\A1(kd1,west))*Dy(north,kd2));
S21=S21+kron(sparse(y1==x2),-(A1(kd1,kd1)\A1(kd1,east))*Dy(south,kd2));
S21=S21+kron(sparse(y1==x1),-(A1(kd1,kd1)\A1(kd1,west))*Dy(south,kd2));

id1=1:(n-2)*size(adjx,1);
id2=(n-2)*size(adjx,1)+(1:(m-2)*size(adjy,1));

sz=(n-2)*size(adjx,1)+(m-2)*size(adjy,1);
S=sparse(sz,sz);
S(id1,id1)=S11;
S(id1,id2)=S12;
S(id2,id1)=S21;
S(id2,id2)=S22;

H=sparse(sz,sz);
H(id1,id1)=H11;
H(id2,id2)=H22;
end