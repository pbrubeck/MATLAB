function [] = HelmTiles( N, k )
e=9;
% Right-Left-Top-Bottom
net=[e 1 e 8; 1 2 e e; 2 e e 3; e e 3 4; 5 e 4 e; 6 5 e e; e 6 7 e; e e 8 7];
RL=[1 2 1 e 2; 2 3 2 1 e; 6 5 5 6 e; 7 6 6 e 5];
TB=[3 4 3 e 4; 4 5 4 3 e; 8 7 7 8 e; 1 8 8 e 7];

% Schur Decompositon Method for the Helmholtz equation on a square with a
% square removed
[D,x]=chebD(N);
D2=D*D;

% Poisson solver on the square [-1/2,1/2]^2
[V,L]=eig(D2(2:N-1,2:N-1), 'vector'); W=inv(V);
[L1,L2]=ndgrid(L); LL=L1+L2;
function [u]=poissonSquare(F)
    u=V*((W*F*W')./LL)*V';
end


% Imposition of Neumann BCs, matching normal derivates at the interface
% The Schur Complement Method maps Dirichlet BCs to Neumann BCs
function [s, h]=schurComplement(b)
    b=reshape(b, N-2, []);
    b=[b, zeros(N-2,1)];
    w=cell([8,1]);
    for j=1:8
        if(norm(b(:,net(j,:)),'fro')==0)
            w{j}=0;
        else
            w{j}=poissonSquare(b(:,net(j,3:4))*D2(2:N-1,[1,N])'+D2(2:N-1,[1,N])*b(:,net(j,1:2))');
        end
    end
    h=zeros(N-2,8);
    s=zeros(N-2,8);
    for j=1:4
        h(:,TB(j,3))=b(:,TB(j,3:5))*[D(N,N)-D(1,1); D(N,1); -D(1,N)];
        h(:,RL(j,3))=[D(N,N)-D(1,1), D(N,1), -D(1,N)]*b(:,RL(j,3:5))';
        s(:,TB(j,3))=h(:,TB(j,3))-(w{TB(j,1)}*D(N,2:N-1)'-w{TB(j,2)}*D(1,2:N-1)');
        s(:,RL(j,3))=h(:,RL(j,3))'-(D(N,2:N-1)*w{RL(j,1)}-D(1,2:N-1)*w{RL(j,2)});
    end
    s=s(:);
    h=h(:);
end


% Construct the Schur complement matrix
S=eye(8*(N-2));
H=eye(8*(N-2));
for i=1:size(S,2)
    [S(:,i), H(:,i)]=schurComplement(S(:,i));
end

% Poisson solver
function [u]=poissonTiles(F)
    F=reshape(F, N-2, N-2, []);
    v=cell([8,1]);
    for j=1:8
        v{j}=poissonSquare(F(:,:,j));
    end
    
    rhs=zeros(N-2,8);
    for j=1:4
        rhs(:,TB(j,3))=-(v{TB(j,1)}*D(N,2:N-1)'-v{TB(j,2)}*D(1,2:N-1)');
        rhs(:,RL(j,3))=-(D(N,2:N-1)*v{RL(j,1)}-D(1,2:N-1)*v{RL(j,2)});
    end
    rhs=rhs(:);

    % Solve for boundary nodes
    b=reshape(S\rhs, N-2, []);
    b=[b, zeros(N-2,1)];
    
    % Solve for interior nodes with the given BCs
    u=zeros(size(F));
    for j=1:8
        u(:,:,j)=poissonSquare(F(:,:,j)-b(:,net(j,3:4))*D2(2:N-1,[1,N])'-D2(2:N-1,[1,N])*b(:,net(j,1:2))');
    end
    u=u(:);
end


% Compute eigenmodes using Arnoldi iteration
[U,lam]=eigs(@poissonTiles, 8*(N-2)^2, k, 'sm');
[lam,id]=sort(diag(lam),'descend');
u=reshape(U(:,id(k)), N-2, N-2, []);
disp(lam/pi^2);

% Retrieve boundary nodes, since they were lost in the eigenmode computation
rhs=zeros(N-2,8);
for i=1:4
    rhs(:,TB(i,3))=-(u(:,:,TB(i,1))*D(N,2:N-1)'-u(:,:,TB(i,2))*D(1,2:N-1)');
    rhs(:,RL(i,3))=-(D(N,2:N-1)*u(:,:,RL(i,1))-D(1,2:N-1)*u(:,:,RL(i,2)));
end
rhs=rhs(:);
b=reshape(H\rhs, N-2, []);



% Plot solution
[xx, yy]=ndgrid([x+2; x(2:end-1); x-2]);
uu=zeros(3*(N-1)+1);
uu(2:end-1,2:end-1)=[u(:,:,1) b(:,8) u(:,:,8) b(:,7) u(:,:,7); 
                     b(:,1)'        zeros(1,N)       b(:,6)';
                     u(:,:,2)      zeros(N-2,N)      u(:,:,6); 
                     b(:,2)'        zeros(1,N)       b(:,5)';
                     u(:,:,3) b(:,3) u(:,:,4) b(:,4) u(:,:,5)];
uu=uu/max(uu(:));

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