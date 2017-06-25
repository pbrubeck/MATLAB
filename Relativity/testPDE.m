function [] = testPDE(m, n)
% Solves a second-order, self-adjoint differential eqation
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
b=[0,0;0,0];
B1=diag(a(1,:))*E1([1,m],:)+diag(b(1,:))*Dx([1,m],:);
B2=diag(a(2,:))*E2([1,n],:)+diag(b(2,:))*Dy([1,n],:);


% Equation coefficients
A11=ones(m,n);
A21=0*ones(m,n)/2.*cos(pi*(3*xx+3*yy));
A12=0*ones(m,n)/2.*cos(pi*(3*xx+3*yy));
A22=ones(m,n);
A0=10*cos(pi*(3*xx+3*yy));

if(any((A12+A21).^2>=4*A11.*A22))
    % This costed a day of work
    error('Equation is not elliptic');
end

% Differential operator
ell=@(uu,ux,uy) Dx*(A11.*ux+A12.*uy)+(A21.*ux+A22.*uy)*Dy'+A0.*uu;
opA=@(uu) ell(uu, Dx*uu, uu*Dy');


% Right hand sides
xb=[x+1i, x-1i];
yb=[1+1i*y; -1+1i*y];
b1=diag(a(1,:))*f(yb)+diag(b(1,:))*fx(yb);
b2=f(xb)*diag(a(2,:))+fy(xb)*diag(b(2,:));
F=opA(f(xx+1i*yy));


% Preconditioner
function b=ophat(A,B,C,D,E,x,tflag)
    X=reshape(x,floor(sqrt(numel(x))),[]);
    if strcmp(tflag,'transp')    
       b=reshape(D*diag(C*diag(B*X*E))*A,[],1);
    else
       b=reshape(B'*diag(C'*diag(D'*X*A'))*E',[],1);
    end
end
function op2=setop(op1,A,B,C,D,E)
    op2=@(x,tflag) op1(x,tflag)+ophat(A,B,C,D,E,x,tflag);
end
ahat=@(x,tflag) ophat(1,1,A0,1,1,x,tflag);
ahat=setop(ahat,Dx,1,A11,Dx,1);
ahat=setop(ahat,1,Dy,A12,Dx,1);
ahat=setop(ahat,Dx,1,A21,1,Dy);
ahat=setop(ahat,1,Dy,A22,1,Dy);

[U,S,V]=svds(ahat, [n*n, m*m], 2);
s=diag(S);
PX1=sqrt(s(1))*reshape(V(:,1),[m,m]);
PX2=sqrt(s(2))*reshape(V(:,2),[m,m]);
PY1=sqrt(s(1))*reshape(U(:,1),[n,n])';
PY2=sqrt(s(2))*reshape(U(:,2),[n,n])';

figure(2);
subplot(2,2,1); imagesc(PX1/sqrt(s(1))); title('PX1'); colormap(gray(32)); colorbar;
subplot(2,2,2); imagesc(PX2/sqrt(s(2))); title('PX2'); colormap(gray(32)); colorbar;
subplot(2,2,3); imagesc(PY1/sqrt(s(1))); title('PY1'); colormap(gray(32)); colorbar;
subplot(2,2,4); imagesc(PY2/sqrt(s(2))); title('PY2'); colormap(gray(32)); colorbar;
drawnow;

% Solver
[gf,ps,kd,gb]=elliptic(PX1,PY1,PX2,PY2,B1,B2,[1,m],[1,n]);

afun=@(uu)  kd(opA(gb(uu)));
pfun=@(rhs) kd(gf(rhs));
ub=ps(b1,b2);
rhs=kd(F-opA(ub));
uu=kd(ub+gf(rhs));
tol=2*eps;
maxit=50;

[uu,~,res,its]=bicgstab(afun,rhs,tol,maxit,pfun,[],uu);
uu=gb(uu)+ub;

display(its);
display(res);

% Interpolation
xq=linspace(-1,1,m);
yq=linspace(-1,1,n);
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