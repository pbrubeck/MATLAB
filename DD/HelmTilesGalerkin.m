function [lam] = HelmTilesGalerkin( m, k )
% Solves the Helmholtz equation using Legedre collocation - weak Galerkin
% spectral method and the Schur complement method for the domain
% decomposition.
n=m;

% UIUC block-I
adjx=[3 4; 4 5; 6 7; 7 8];
adjy=[4 2; 2 1; 1 7];

% L-shaped membrane
adjx=[2 1];
adjy=[3 1];

% Topology
[topo,net,RL,TB]=ddtopo(adjx,adjy);
pos=ddpatches(topo);
x0=real(pos);
y0=imag(pos);

% Degrees of freedom
rd1=[1,m];
rd2=[1,n];
kd1=2:m-1;
kd2=2:n-1;

% Schur complement
[S,H,gf,K1,K2,M1,M2,x,y]=ddschurGalerkin(adjx,adjy,m,n,ones(2,2),zeros(2,2));
[Lschur, Uschur, pschur]=lu(S,'vector');

figure(2);
imagesc(log(abs(S)));
colormap(gray(256)); colorbar; axis square;
drawnow;


% Poisson solver
function [u]=poissonTiles(F)
    F=reshape(F, m, n, []);
    v=cell([size(net,1),1]);
    for j=1:size(net,1)
        v{j}=gf(F(:,:,j),zeros(2,n-2),zeros(m-2,2));
    end
    
    rhs=zeros(m-2, size(RL,1)+size(TB,1));
    for j=1:size(RL,1)
        rhs(:,RL(j,1))=M1(rd1(2),:)*F(:, :, adjx(j,1))*M2(kd2,:)'/2 + M1(rd1(1),:)*F(:, :, adjx(j,2))*M2(kd2,:)'/2 + ...
                       -(M1(rd1(2),:)*(v{adjx(j,1)})*K2(kd2,:)'+K1(rd1(2),:)*(v{adjx(j,1)})*M2(kd2,:)' + ...
                         M1(rd1(1),:)*(v{adjx(j,2)})*K2(kd2,:)'+K1(rd1(1),:)*(v{adjx(j,2)})*M2(kd2,:)');
    end
    for j=1:size(TB,1)
        rhs(:,TB(j,1))=M1(kd1,:)*F(:, :, adjy(j,1))*M2(rd2(2),:)'/2 + M1(kd1,:)*F(:, :, adjy(j,2))*M2(rd2(1),:)'/2 + ...
                       -(K1(kd1,:)*(v{adjy(j,1)})*M2(rd2(2),:)'+M1(kd1,:)*(v{adjy(j,1)})*K2(rd2(2),:)'+...
                         K1(kd1,:)*(v{adjy(j,2)})*M2(rd2(1),:)'+M1(kd1,:)*(v{adjy(j,2)})*K2(rd2(1),:)');
    end
    rhs=rhs(:);
    
    % Solve for boundary nodes
    b=Uschur\(Lschur\rhs(pschur));
    b=reshape(b, m-2, []);
    b=[b, zeros(m-2,1)];
    
    % Solve for interior nodes with the given BCs
    u=zeros(size(F));
    for j=1:size(net,1)
        u(:,:,j)=gf(F(:,:,j), b(:,net(j,1:2))', b(:,net(j,3:4)));
    end
    u=u(:);   
end

% Compute eigenmodes using Arnoldi iteration
[U,lam]=eigs(@poissonTiles, size(net,1)*(m*n), k, 'sm');
[lam,id]=sort(real(diag(lam)),'ascend');
U=U(:,id);

um=reshape(U, m, n, [], k);
uuu=permute(um,[1,2,4,3]);

[xx,yy]=ndgrid(x,y);

figure(1); zoom off; pan off; rotate3d off;
for i=1:size(uuu,4)
    modegallery(xx+x0(i),yy+y0(i),uuu(:,:,:,i));
    if i==1, hold on; end;
end
hold off;
colormap(jet(256));  shading interp; camlight; view(2);
xl=xlim(); dx=xl(2)-xl(1);
yl=ylim(); dy=yl(2)-yl(1);
pbaspect([dx,dy,min(dx,dy)]);
xlabel('x'); ylabel('y');
end