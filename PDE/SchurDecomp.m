function [lam] = SchurDecomp(N, k)
% Schur Decompositon Method for the Helmholtz equation on the L-shaped
% membrane

% Poisson solver on the square [0,1]^2
[D,x]=chebD(N); x=x/2; D=2*D; D2=D*D;
[V,L]=eig(D2(2:N-1,2:N-1), 'vector'); W=inv(V);
[L1,L2]=ndgrid(L); LL=L1+L2;
function [u]=poissonSquare(F)
    u=V*((W*F*W')./LL)*V';
end


% Imposition of Neumann BCs, matching normal derivates at the interface
% The Schur Complement Method maps Dirichlet BCs to Neumann BCs
function [rhs]=schurComplement(b)
    b=reshape(b, N-2, []);
    w1=poissonSquare(D2(2:N-1,N)*b(:,1)');
    w2=poissonSquare(D2(2:N-1,1)*b(:,1)'+b(:,2)*D2(2:N-1,N)');
    w3=poissonSquare(b(:,2)*D2(2:N-1,1)');
    r1=(D(N,N)-D(1,1))*b(:,1)'-(D(N,2:N-1)*w1-D(1,2:N-1)*w2);
    r2=(D(N,N)-D(1,1))*b(:,2)-(w2*D(N,2:N-1)'-w3*D(1,2:N-1)');
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
    
    rhs1=-(D(N,2:N-1)*v1-D(1,2:N-1)*v2);
    rhs2=-(v2*D(N,2:N-1)'-v3*D(1,2:N-1)');
    rhs=[rhs1(:); rhs2(:)];

    % Solve for boundary nodes
    b=reshape(S\rhs, N-2, []);

    % Solve for interior nodes with the given BCs
    u1=poissonSquare(F(:,:,1)-D2(2:N-1,N)*b(:,1)');
    u2=poissonSquare(F(:,:,2)-D2(2:N-1,1)*b(:,1)'-b(:,2)*D2(2:N-1,N)');
    u3=poissonSquare(F(:,:,3)-b(:,2)*D2(2:N-1,1)');
    u=[u1(:); u2(:); u3(:)];
end


% Compute eigenmodes using Arnoldi iteration
[U,lam]=eigs(@poissonL, 3*(N-2)^2, k, 'sm');
[lam,id]=sort(diag(lam),'descend');
u=reshape(U(:,id(k)), N-2, N-2, []);

% Retrieve boundary nodes, since they were lost in the eigenmode computation
b0=zeros(N-2,1);
b1=-(D(N,2:N-1)*u(:,:,1)-D(1,2:N-1)*u(:,:,2))/(D(N,N)-D(1,1));
b2=-(u(:,:,2)*D(N,2:N-1)'-u(:,:,3)*D(1,2:N-1)')/(D(N,N)-D(1,1));

umax=max(u(:));
u=u/umax;
b1=b1/umax;
b2=b2/umax;

% Plot solution
[xx, yy]=ndgrid(x);

figure(1);
surf(xx+1/2, yy+1/2, [0 b0' 0; b0 u(:,:,1) b0; 0 b1  0]); hold on;
surf(xx-1/2, yy+1/2, [0 b1  0; b0 u(:,:,2) b2; 0 b0' 0]); 
surf(xx-1/2, yy-1/2, [0 b0' 0; b2 u(:,:,3) b0; 0 b0' 0]); hold off;

colormap(jet(256)); view(2); shading interp;
title(sprintf('\\lambda_{%d}=%.8f', k, lam(k)));
xlabel('x'); ylabel('y'); axis square; axis manual;
end