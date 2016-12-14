function [] = schurDecomp( N, k )
% Schur Decompositon Method for the Helmholtz equation on the L-shaped
% membrane

[D,x]=chebD(N); x=x/2; D=2*D;
D2=D*D;

[V,L]=eig(D2(2:N-1,2:N-1), 'vector');
W=inv(V);
[Lx,Ly]=ndgrid(L);
LL=Lx+Ly;


% Poisson solver on the square [-1/2,1/2]^2
function [u]=poissonSquare(F)
    u=V*((W*F*W')./LL)*V';
end


% Imposition of Neumann BCs, matching normal derivates at the interface
% The Schur Complement Method maps Dirichlet BCs to Neumann BCs
function [rhs]=schurComplement(b)
    b=reshape(b,N-2,[]);
    w1=poissonSquare(b(:,1)*D2(2:N-1,1)'+D2(2:N-1,1)*b(:,2)');
    w2=poissonSquare(b(:,1)*D2(2:N-1,N)');
    w3=poissonSquare(D2(2:N-1,N)*b(:,2)');
    r1=(D(end)-D(1))*b(:,1)-(w2*D(N,2:N-1)'-w1*D(1,2:N-1)');
    r2=(D(end)-D(1))*b(:,2)'-(D(N,2:N-1)*w3-D(1,2:N-1)*w1);
    rhs=[r1(:); r2(:)];
end


% Construct the Schur complement matrix
S=eye(2*(N-2));
for i=1:size(S,2)
    S(:,i)=schurComplement(S(:,i));
end

% Poisson solver on the L-shaped membrane
function [u]=poissonL(F)
    F=reshape(F, N-2, N-2, []);
    v1=poissonSquare(F(:,:,1));
    v2=poissonSquare(F(:,:,2));
    v3=poissonSquare(F(:,:,3));

    rhs1=-(v2*D(N,2:N-1)'-v1*D(1,2:N-1)');
    rhs2=-(D(N,2:N-1)*v3-D(1,2:N-1)*v1);
    rhs=[rhs1(:); rhs2(:)];

    % Solve for boundary nodes
    b=S\rhs;
    b1=b(1:N-2); b2=b(N-1:2*(N-2))';

    % Solve for interior nodes with the given BCs
    u1=poissonSquare(F(:,:,1)-D2(2:N-1,1)*b2-b1*D2(2:N-1,1)');
    u2=poissonSquare(F(:,:,2)-b1*D2(2:N-1,N)');
    u3=poissonSquare(F(:,:,3)-D2(2:N-1,N)*b2);
    u=[u1(:); u2(:); u3(:)];
end


% Compute eigenmodes using Arnoldi iteration
[U,lam]=eigs(@poissonL, 3*(N-2)^2, k, 'sm');
u=reshape(U(:,k), N-2, N-2, []);
lam=diag(lam);

% Retrieve boundary nodes, since they were lost in the eigenmode computation
b1=-(u(:,:,2)*D(N,2:N-1)'-u(:,:,1)*D(1,2:N-1)')/(D(end)-D(1));
b2=-(D(N,2:N-1)*u(:,:,3)-D(1,2:N-1)*u(:,:,1))/(D(end)-D(1));

% Plot solution
uu=zeros(2*N-1);
uu(2:end-1,2:end-1)=[zeros(N-2, N-1), u(:,:,3); zeros(1,N-1), b2; u(:,:,2), b1, u(:,:,1)];
[xx, yy]=ndgrid([x+1/2; x(2:end)-1/2]);


figure(1);
surfl(xx, yy, uu, 'light');
title(sprintf('\\lambda_{%d}=%f', k, lam(k)));
zrange=max(uu(:))-min(uu(:));
xrange=max(xx(:))-min(xx(:));
yrange=max(yy(:))-min(yy(:));
daspect([1 1 2*zrange/hypot(xrange,yrange)]);
xlabel('x'); ylabel('y');
colormap(jet(256)); view(2); shading interp;
end