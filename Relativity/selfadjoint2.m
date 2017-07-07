function [] = selfadjoint2(m, n)
% Solves a second-order, self-adjoint differential equation
if nargin<2
    n=m;
end

% Exact solution
f=@(z) real(sin(z.^4));
fx=@(z) real(4*z.^3.*cos(z.^4));
fy=@(z) real(1i*4*z.^3.*cos(z.^4));

% Differential operators
E1=eye(m);
E2=eye(n);
[Dx,x]=chebD(m);
[Dy,y]=chebD(n); y=y';
[xx,yy]=ndgrid(x,y);


% Boundary condtions
a=[1,1;1,1];
b=0*[1,-1;1,-1];
B1=diag(a(1,:))*E1([1,m],:)+diag(b(1,:))*Dx([1,m],:);
B2=diag(a(2,:))*E2([1,n],:)+diag(b(2,:))*Dy([1,n],:);


% Equation coefficients
[G11,G21,G12,G22,K11,K21,K12,K22,vol]=mymanifold(x,y,'inv');
%[map,vol,G11,G12,G22]=mapquad([1+1i;2-1i;-0.9+0.8i;-1.3-0.7i],'inv');

A11=vol.*G11;
A21=vol.*G12;
A12=vol.*G12;
A22=vol.*G22;
A0=vol.*zeros(m,n);

figure(2);
subplot(2,2,1); surf(xx,yy,A11); title('A11'); shading interp; colormap(jet(256)); colorbar; view(2);
subplot(2,2,2); surf(xx,yy,A12); title('A12'); shading interp; colormap(jet(256)); colorbar; view(2);
subplot(2,2,3); surf(xx,yy,A21); title('A21'); shading interp; colormap(jet(256)); colorbar; view(2);
subplot(2,2,4); surf(xx,yy,A22); title('A22'); shading interp; colormap(jet(256)); colorbar; view(2);

if(any((A12+A21).^2>=4*A11.*A22))
    % This costed a day of work
    error('Equation is not elliptic');
end

% Self-adjoint differential operator
ell=@(uu,ux,uy) Dx*(A11.*ux+A12.*uy)+(A21.*ux+A22.*uy)*Dy'+A0.*uu;
opA=@(uu) ell(uu, Dx*uu, uu*Dy');

% Right hand sides
zz=xx+1i*yy;
z1=zz([1,end],:);
z2=zz(:,[1,end]);
b1=diag(a(1,:))*f(z1)+diag(b(1,:))*fx(z1);
b2=f(z2)*diag(a(2,:))+fy(z2)*diag(b(2,:));
F=opA(f(xx+1i*yy));


% Preconditioner
function b=ophat(A,B,C,D,E,x,tflag)
    X=reshape(x,floor(sqrt(numel(x))),[]);
    if strcmp(tflag,'notransp')
        b=reshape(E*diag(diag(A*X.'*D).'*C)*B,[],1);
    else
        b=reshape(D*diag(C*diag(B*X.'*E))*A,[],1);
    end
end
function op2=setop(op1,A,B,C,D,E)
    op2=@(x,tflag) op1(x,tflag)+ophat(A,B,C,D,E,x,tflag);
end

% Rank-1 approximation seems to be optimal for A11 and A22
function B=lowrank(A,r)
    [U0,S0,V0]=svds(A,r);
    B=U0*S0*V0';
end

ahat=@(x,tflag) ophat(1,1,A0,1,1,x,tflag);
ahat=setop(ahat,Dx,1,lowrank(A11,1),Dx,1);
ahat=setop(ahat,1,Dy,A12,Dx,1);
ahat=setop(ahat,Dx,1,A21,1,Dy);
ahat=setop(ahat,1,Dy,lowrank(A22,1),1,Dy);

[U,S,V]=svds(ahat, [n*n, m*m], 2);
s=diag(S);
PX1=sqrt(s(1))*reshape(V(:,1),[m,m]);
PX2=sqrt(s(2))*reshape(V(:,2),[m,m]);
PY1=sqrt(s(1))*reshape(U(:,1),[n,n]);
PY2=sqrt(s(2))*reshape(U(:,2),[n,n]);

figure(3);
subplot(2,2,1); imagesc(PX1/sqrt(s(1))); title('PX1'); colormap(gray(32)); colorbar;
subplot(2,2,2); imagesc(PX2/sqrt(s(2))); title('PX2'); colormap(gray(32)); colorbar;
subplot(2,2,3); imagesc(PY1/sqrt(s(1))); title('PY1'); colormap(gray(32)); colorbar;
subplot(2,2,4); imagesc(PY2/sqrt(s(2))); title('PY2'); colormap(gray(32)); colorbar;

% Solver
[gf,ps,kd,gb]=elliptic(PX1,PY1,PX2,PY2,B1,B2,[1,m],[1,n]);

afun=@(uu)  kd(opA(gb(uu)));
pfun=@(rhs) kd(gf(rhs));
ub=ps(b1,b2);
rhs=kd(F-opA(ub));
uu=kd(ub+gf(rhs));
tol=2*eps;
maxit=25;

[uu,~,res,its]=bicgstab(afun,rhs,tol,maxit,pfun,[],uu);
uu=gb(uu)+ub;

display(its);
display(res);

% Eigenvalue solver
solver=@(rhs) bicgstab(afun,rhs,tol,maxit,pfun,[],pfun(rhs));
[uu,L]=eigs(solver, numel(rhs), 25, 'sm');
[L,id]=sort(diag(real(L)),'descend');
uu=gb(uu(:,id(end)));
display(L);

% Interpolation
xq=linspace(-1,1,1024);
yq=linspace(-1,1,1024);
[xxx,yyy]=ndgrid(xq, yq);
uuu=interpcheb(interpcheb(uu,xq,1),yq,2);

err=norm(uu-f(xx+1i*yy),'inf');
errint=norm(uuu-f(xxx+1i*yyy),'inf');
display(err);
display(errint);

% Plot
figure(1);
surf(xxx,yyy,uuu);
colormap(jet(256));
shading interp;
camlight;
end