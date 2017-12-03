function [lam] = HelmTiles( N, k )
% Schur Decompositon Method for the Helmholtz equation on a domain formed
% by squares

adjx=[2 3; 3 4; 5 6; 6 7];
adjy=[3 1; 1 6];
[topo,net,RL,TB]=ddtopo(adjx,adjy);
pos=ddpatches(topo);
x0=real(pos);
y0=imag(pos);

% Poisson solver on the square [-1,1]^2
[D,x]=chebD(N);
D2=D*D;
[V,L]=eig(D2(2:N-1,2:N-1), 'vector'); W=inv(V);
[L1,L2]=ndgrid(L); LL=L1+L2;
poissonSquare=@(F) V*((W*F*W')./LL)*V';

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
    for j=1:size(RL,1)
        h(:,RL(j,1))=[D(N,N)-D(1,1), D(N,1), -D(1,N)]*b(:,RL(j,:))';
        s(:,RL(j,1))=h(:,RL(j,1))'-(D(N,2:N-1)*w{adjx(j,1)}-D(1,2:N-1)*w{adjx(j,2)});
    end
    for j=1:size(TB,1)
        h(:,TB(j,1))=b(:,TB(j,:))*[D(N,N)-D(1,1); D(N,1); -D(1,N)];
        s(:,TB(j,1))=h(:,TB(j,1))-(w{adjy(j,1)}*D(N,2:N-1)'-w{adjy(j,2)}*D(1,2:N-1)');
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
    for j=1:size(RL,1)
        rhs(:,RL(j,1))=-(D(N,2:N-1)*v{adjx(j,1)}-D(1,2:N-1)*v{adjx(j,2)});
    end
    for j=1:size(TB,1)
        rhs(:,TB(j,1))=-(v{adjy(j,1)}*D(N,2:N-1)'-v{adjy(j,2)}*D(1,2:N-1)');
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
U=U(:,id);

um=reshape(U, N-2, N-2, [], k);
uuu=zeros(N,N,k,size(um,3));
[xx, yy]=ndgrid(x);
for mode=1:k
    % Retrieve boundary nodes, since they were lost in the eigenmode computation
    rhs=zeros(N-2,size(RL,1)+size(TB,1));
    for i=1:size(RL,1)
        rhs(:,RL(i,1))=-(D(N,2:N-1)*um(:,:,adjx(i,1),mode)-D(1,2:N-1)*um(:,:,adjx(i,2),mode));
    end
    for i=1:size(TB,1)
        rhs(:,TB(i,1))=-(um(:,:,adjy(i,1),mode)*D(N,2:N-1)'-um(:,:,adjy(i,2),mode)*D(1,2:N-1)');
    end
    b=[reshape(H\rhs(:), N-2, []), zeros(N-2,1)];
    % Assemble mode
    for i=1:size(um,3)
        uuu(2:N-1,2:N-1,mode,i)=um(:,:,i,mode);
        uuu([1,N],2:N-1,mode,i)=b(:,net(i,1:2))';
        uuu(2:N-1,[1,N],mode,i)=b(:,net(i,3:4));
    end
end

figure(1);
for i=1:size(uuu,4)
    modegallery(xx+x0(i),yy+y0(i),uuu(:,:,:,i));
    if i==1, hold on; end;
end
hold off;
colormap(jet(256)); camlight; view(2); shading interp; 
title(sprintf('\\lambda_{%d}=%.8f', k, lam(k)));
xlabel('x'); ylabel('y');
end