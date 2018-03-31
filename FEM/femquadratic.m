function [err,h] = femquadratic(gmshfile,n)
% Canonical domain
z0=[0; 1; 1i; 0.5; 0.5+0.5i; 0.5i];
x0=real(z0);
y0=imag(z0);
P0=[ones(size(x0)), 2*x0, 2*y0, x0.^2, 2*x0.*y0, y0.^2];
R=full(sparse(1:9,[1,2,3,2,4,5,3,5,6],1,9,6));

% Permutation operator (maps [a,b;c,d] to [d,-b;-c,a])
Perm=full(sparse(1:4,[4,2,3,1],[1,-1,-1,1],4,4));

% Quadrature in canonical domain
[zq,wq]=trigauss(n); 

Vq=[ones(size(zq)), real(zq), imag(zq)];

% Interpolation to quadrature points, evaluate quadratic
H=R/P0;
E=zeros(length(wq), 6);
for ej=1:6
    E(:,ej)=sum((Vq*reshape(H(:,ej),[3,3])).*Vq,2);
end

% Derivative of cardinal basis functions, evaluate linear
Dx=full(sparse(1:3,[2,4,5],2,3,6));
Dy=full(sparse(1:3,[3,5,6],2,3,6));
Dx=Vq*(Dx/P0);
Dy=Vq*(Dy/P0);


[tri, vert, bnd] = loadgmsh(gmshfile);
N=size(vert,1);
T=size(tri,1);

h=0;
m=6;
ei=zeros(m*m*T,1);
ej=zeros(m*m*T,1);
Mass=zeros(m*m*T,1);
Stiff=zeros(m*m*T,1);
for e=1:T    
    nodes=tri(e,:);
    Pe=[ones(3,1), real(vert(nodes(1:3))), imag(vert(nodes(1:3)))];
    area=abs(det(Pe))/2;
    h=max(h, sqrt(2*area));
    
    % Element vertices
    ze=vert(nodes);

    % Element mapping
    He=reshape(R*(P0\ze), 3,3);
    Je=2*He(:,2:3);
    % w=sum((Vq*He).*Vq,2);

    % Element metric
    Ge=real((kron(Vq,ones(1,3)).*kron(ones(1,3),Vq))*kron(conj(Je),Je));
    jac=sqrt(Ge(:,1).*Ge(:,4)-Ge(:,2).*Ge(:,3));
    % Invert metric
    Ge = (Ge*Perm)./jac;

    % Element mass and stiffness
    Me=E'*((wq.*jac).*E);
    Ke=Dx'*((wq.*Ge(:,1)).*Dx+(wq.*Ge(:,2)).*Dy) + ...
       Dy'*((wq.*Ge(:,3)).*Dx+(wq.*Ge(:,4)).*Dy);

    idx=(m*m)*(e-1)+1:(m*m)*e;
    [ni,nj]=ndgrid(nodes);
    ei(idx)=ni(:);
    ej(idx)=nj(:);
    Mass(idx)=Me(:);  % add Me to global Mass
    Stiff(idx)=Ke(:); % add Ke to global Stiff
end

natural=find(((real(vert)==0) | (imag(vert)==0)) ...
    & ~(vert==1) & ~(vert==2) & ~(vert==1i) & ~(vert==2i));
bnd=setdiff(bnd,natural);

er=1:N;
er(bnd)=[];
Rassembly=sparse(er,1:numel(er),1,N,numel(er));

% Assembly
Mass =sparse(ei, ej, Mass, N, N);
Stiff=sparse(ei, ej, Stiff, N, N);

% Poisson right hand side
f=ones(size(vert));

% Solving
rhs=Rassembly'*(Mass*f(:));
u=Rassembly*((Rassembly'*Stiff*Rassembly)\rhs);

% Exact solution
r=abs(vert);
uex = (1-r.^2)/4 + 3/4*log(r)/log(2);
err = norm(u-uex,'inf')/norm(uex,'inf');

% Ploting
tri2=delaunay(real(vert),imag(vert));
kd=zeros(size(tri2,1));
for e=1:size(tri2,1)
    nodes=tri2(e,:);
    s=max(abs([0,1,-1;1,0,-1;1,-1,0]*vert(nodes)));
    if(s<2*h) 
        kd(e)=e; 
    end
end
tri2=tri2(kd(kd>0),:);

figure(1);
trisurf(tri2,real(vert),imag(vert),u,u);
colormap(jet(256)); colorbar('TickLabelInterpreter','latex');
axis square; view(2);
alpha(0.8);

set(gcf,'defaultTextInterpreter','latex');
set(gca,'TickLabelInterpreter','latex','fontsize',14);
end





function [tri, vert, bnd] = loadmesh(vfile, efile)
fileID = fopen(vfile,'r');
vert = fscanf(fileID, '%f %f', [2 Inf])';
fclose(fileID);

fileID = fopen(efile,'r');
tri = fscanf(fileID, '%d %d %d', [3 Inf])';
fclose(fileID);

tri=tri+1;
vert=vert*[1;1i];


sides=[1,2;2,3;1,3];
p=0;
qvert=zeros([3*size(tri,1),1]); % Overkill
qtri=zeros(size(tri));
% First pass, element interfaces
for i=1:size(tri,1)
    j=i+1;
    while(j<=size(tri,1) && ~all(qtri(i,:)))
        for s1=find(qtri(i,:)==0)
            pair1 = tri(i, sides(s1,:));
            for s2=1:size(sides,1)
                pair2 = tri(j, sides(s2,:));
                if(all(pair1==pair2) || all(pair1==pair2([2,1])))
                    p=p+1;
                    mid=p+size(vert,1);
                    qtri(i,s1)=mid;
                    qtri(j,s2)=mid;
                    qvert(p)=mean(vert(pair1));
                end
            end
        end
        j=j+1;
    end
end

bnd=zeros(size(vert,1)+size(qvert,1),1);
% Second pass, boundaries
for i=1:size(tri,1)
    for s1=find(qtri(i,:)==0)
        pair1=tri(i, sides(s1,:));
        p=p+1;
        mid=p+size(vert,1);
        qtri(i,s1)=mid;
        qvert(p)=mean(vert(pair1));
        bnd([pair1, mid])=1;
    end
end
qvert=qvert(1:p);
vert=[vert; qvert];
tri=[tri, qtri];
bnd=find(bnd);
end