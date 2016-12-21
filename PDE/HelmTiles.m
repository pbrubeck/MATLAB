function [lam] = HelmTiles( N, k )
% Schur Decompositon Method for the Helmholtz equation on a domain formed
% by squares

% Inner square cavity
x0=[2;0;-2;-2;-2;0;2;2];
y0=[2;2;2;0;-2;-2;-2;0];
e=9;
% Right-Left-Top-Bottom
net=[e 1 e 8; 1 2 e e; 2 e e 3; e e 3 4; 5 e 4 e; 6 5 e e; e 6 7 e; e e 8 7];
% Id-NextB-LastB-NextU-LastU
RL=[1 e 2 1 2; 2 1 e 2 3; 5 6 e 6 5; 6 e 5 7 6];
TB=[3 e 4 3 4; 4 3 e 4 5; 7 8 e 8 7; 8 e 7 1 8];

% L-shaped membrane
x0=[1;-1;-1];
y0=[1;1;-1];
e=3;
net=[e 1 e e; 1 e e 2; e e 2 e];
RL=[1 e e 1 2];
TB=[2 e e 2 3];

% Poisson solver on the square [-1,1]^2
[D,x]=chebD(N);
D2=D*D;
[V,L]=eig(D2(2:N-1,2:N-1), 'vector'); W=inv(V);
[L1,L2]=ndgrid(L); LL=L1+L2;
function [u]=poissonSquare(F)
    u=V*((W*F*W')./LL)*V';
end


% Imposition of Neumann BCs, matching normal derivates at the interface
% The Schur Complement Method maps Dirichlet BCs to Neumann BCs
function [s, h]=schurComplement(b)
    b=reshape(b, N-2, []);
    h=zeros(size(b));
    s=zeros(size(b));
    b=[b, zeros(N-2,1)];
    w=cell([size(net,1),1]);
    for j=1:size(net,1)
        if(any(any(b(:,net(j,:)))))
            w{j}=poissonSquare(D2(2:N-1,[1,N])*b(:,net(j,1:2))'+b(:,net(j,3:4))*D2(2:N-1,[1,N])');
        else
            w{j}=0;
        end
    end
    for j=1:size(TB,1)
        h(:,TB(j,1))=b(:,TB(j,1:3))*[D(N,N)-D(1,1); D(N,1); -D(1,N)];
        s(:,TB(j,1))=h(:,TB(j,1))-(w{TB(j,4)}*D(N,2:N-1)'-w{TB(j,5)}*D(1,2:N-1)');
    end
    for j=1:size(RL,1)
        h(:,RL(j,1))=[D(N,N)-D(1,1), D(N,1), -D(1,N)]*b(:,RL(j,1:3))';
        s(:,RL(j,1))=h(:,RL(j,1))'-(D(N,2:N-1)*w{RL(j,4)}-D(1,2:N-1)*w{RL(j,5)});
    end
    s=s(:);
    h=h(:);
end


% Construct the Schur complement matrix
S=eye((size(RL,1)+size(TB,1))*(N-2));
H=eye((size(RL,1)+size(TB,1))*(N-2));
for i=1:size(S,2)
    [S(:,i), H(:,i)]=schurComplement(S(:,i));
end

% Poisson solver
function [u]=poissonTiles(F)
    F=reshape(F, N-2, N-2, []);
    v=cell([size(net,1),1]);
    for j=1:size(net,1)
        v{j}=poissonSquare(F(:,:,j));
    end
    
    rhs=zeros(N-2, size(RL,1)+size(TB,1));
    for j=1:size(TB,1)
        rhs(:,TB(j,1))=-(v{TB(j,4)}*D(N,2:N-1)'-v{TB(j,5)}*D(1,2:N-1)');
    end
    for j=1:size(RL,1)
        rhs(:,RL(j,1))=-(D(N,2:N-1)*v{RL(j,4)}-D(1,2:N-1)*v{RL(j,5)});
    end
    rhs=rhs(:);

    % Solve for boundary nodes
    b=reshape(S\rhs, N-2, []);
    b=[b, zeros(N-2,1)];
    
    % Solve for interior nodes with the given BCs
    u=zeros(size(F));
    for j=1:size(net,1)
        u(:,:,j)=poissonSquare(F(:,:,j)-D2(2:N-1,[1,N])*b(:,net(j,1:2))'-b(:,net(j,3:4))*D2(2:N-1,[1,N])');
    end
    u=u(:);
end


% Compute eigenmodes using Arnoldi iteration
[U,lam]=eigs(@poissonTiles, size(net,1)*(N-2)^2, k, 'sm');
[lam,id]=sort(real(diag(lam)),'descend');
u=reshape(U(:,id(k)), N-2, N-2, []);


% Retrieve boundary nodes, since they were lost in the eigenmode computation
rhs=zeros(N-2,size(RL,1)+size(TB,1));
for i=1:size(TB,1)
    rhs(:,TB(i,1))=-(u(:,:,TB(i,4))*D(N,2:N-1)'-u(:,:,TB(i,5))*D(1,2:N-1)');
end
for i=1:size(RL,1)
    rhs(:,RL(i,1))=-(D(N,2:N-1)*u(:,:,RL(i,4))-D(1,2:N-1)*u(:,:,RL(i,5)));
end
rhs=rhs(:);
b=reshape(H\rhs, N-2, []);
b=[b, zeros(N-2,1)];

% Plot solution
[xx, yy]=ndgrid(x);
umax=max(u(:));
uu=zeros(N);
figure(1);
for i=1:size(u,3)
    uu(2:N-1,2:N-1)=u(:,:,i);
    uu([1,N],2:N-1)=b(:,net(i,1:2))';
    uu(2:N-1,[1,N])=b(:,net(i,3:4));
    surf(xx+x0(i), yy+y0(i), uu/umax);
    if i==1, hold on; end;
end
hold off;

colormap(jet(256)); view(2); shading interp;
title(sprintf('\\lambda_{%d}=%.8f', k, lam(k)));
xlabel('x'); ylabel('y'); axis square manual;
end