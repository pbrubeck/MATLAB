function [] = ngon(n, N, k)
% Dirichlet eigenmodes of the laplacian on regular n-gon mapped to the
% unit disk
p=polygon(exp(2*pi*1i*(0:n-1)/n));
f=diskmap(p);
f=center(f,0);

% Set differential operators
N([1,2])=N;
[R2,Drr,Dtt,r,th]=chebLapPol(N(1),N(2));
R2=R2(2:end,2:end);

% Jacobian determinant
zz=r*exp(1i*th);
J=abs(evaldiff(f,zz(2:end,:))).^2;

% Compute the k-th eigenmode
[V1, w1]=eig(Drr(2:end,2:end),'vector'); U1=inv(V1);
[V2, w2]=eig(Dtt,'vector'); U2=inv(V2);
[W1,W2]=ndgrid(w1,w2); W=W1+W2;
OP=@(x) reshape(greenF(V1,V2,U1,U2,W,R2*J,reshape(x, size(J))), size(x));
[V,lam]=eigs(OP, numel(J), k, 'sm');
lam=diag(lam);
[lam,idx]=sort(lam,'descend');
phi=zeros(N);
phi(2:end,:)=reshape(V(:,idx(k)), size(J));
phi=phi/max(abs(phi(:)));
disp(lam);

% Conformal mapping ww=f(zz)
b1=f(zz(1,:));
BC=Drr(:,1)*b1;
ww=zeros(N); ww(1,:)=b1;
ww(2:end,:)=greenF(V1,V2,U1,U2,W,R2*ones(size(J)),-BC(2:end,:));
uu=real(ww); vv=imag(ww);

phi=[phi(1:N,end), phi(1:N,:)];
uu=[uu(1:N,end), uu(1:N,:)];
vv=[vv(1:N,end), vv(1:N,:)];

% Plot solution
figure(1);
surfl(uu,vv,phi,'light'); 
shading interp; alpha(0.8);
colormap(jet(256));
zrange=max(phi(:))-min(phi(:));
xrange=max(uu(:))-min(uu(:));
yrange=max(vv(:))-min(vv(:));
daspect([1 1 2*zrange/hypot(xrange,yrange)]);
xlim([-1, 1]);
ylim([-1, 1]);
view(0,90);
title(sprintf('\\lambda_{%d} = %.8f', k, lam(k)));
end

function X=greenF(V1, V2, U1, U2, W, J, X)
X=V1*((U1*(J.*X)*U2')./W)*V2';
end