function [] = schurDecomp( N )
% Schur Decompositon Method for the Poisson equation on the L-shaped
% membrane

F1=-ones(N-2);
F2=-ones(N-2);
F3=-ones(N-2);

[D,x]=chebD(N);
D2=D*D;

A=D2(2:N-1,2:N-1);
[V,L]=eig(A,'vector');
W=inv(V);
[Lx,Ly]=ndgrid(L);
LL=Lx+Ly;

function [u]=poissonSquare(F)
    u=V*((W*F*W')./LL)*V';
end

% Imposition of Neumann BCs, matching normal derivates at the interface
% The Schur Complement Method maps Dirichlet BCs to Neumann BCs
function [rhs]=schurComplement(b)
    b1=b(1:N-2);
    b2=b(N-1:2*(N-2))';
    
    T1=poissonSquare(D2(2:N-1,1)*b2+b1*D2(2:N-1,1)');
    T2=poissonSquare(b1*D2(2:N-1,N)');
    T3=poissonSquare(D2(2:N-1,N)*b2);
    
    r1=(D(end)-D(1))*b1-(T2*D(N,2:N-1)'-T1*D(1,2:N-1)');
    r2=(D(end)-D(1))*b2-(D(N,2:N-1)*T3-D(1,2:N-1)*T1);
    rhs=[r1(:); r2(:)];
end

v1=poissonSquare(F1);
v2=poissonSquare(F2);
v3=poissonSquare(F3);

rhs1=-(v2*D(N,2:N-1)'-v1*D(1,2:N-1)');
rhs2=-(D(N,2:N-1)*v3-D(1,2:N-1)*v1);
rhs=[rhs1(:); rhs2(:)];

b=cgs(@schurComplement, rhs, 1e-06, 2*N);
b1=b(1:N-2); b2=b(N-1:2*(N-2))';

u1=v1-T1;
u2=v2-T2;
u3=v3-T3;


% plot solution
uu=zeros(2*N-1);
uu(2:end-1,2:end-1)=[zeros(N-2, N-1), u3; zeros(1,N-1), b2; u2, b1, u1];
[xx, yy]=ndgrid([x+1; x(2:end)-1]);
surf(xx, yy, uu);
zrange=max(uu(:))-min(uu(:));
xrange=max(xx(:))-min(xx(:));
yrange=max(yy(:))-min(yy(:));
daspect([1 1 2*zrange/hypot(xrange,yrange)]);
xlabel('x'); ylabel('y');
colormap(jet(256)); view(2); %shading interp;
end