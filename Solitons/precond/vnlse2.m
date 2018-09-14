function [E] = vnlse2(n)
% Stiffness matrix
[D,z,w]=legD(n);
hx=10;
x=hx*(z);
wx=hx*w;
K1=D'*diag(w/hx)*D;
M1=diag(wx);
B1=diag(wx);

hy=10;
y=hy*(z);
wy=hy*w;
K2=D'*diag(w/hy)*D;
M2=diag(wy);

lam=1/2;
H1=K1+lam*B1;
H2=K2;

[xx,yy]=ndgrid(x,y);
rr=hypot(yy,xx);
tt=atan2(yy,xx);

H=@(psi) reshape( H1*psi*M2' + M1*psi*H2',[],1);

% Nonlinear term
bess=-0*besselj(1,5*rr).^2;
nl=@(psi2) -psi2.^2/2+bess.*psi2;
% Lagrangian
lag=@(psi) (psi(:)'*H(psi) + wx'*nl(psi.*conj(psi))*wy)/2;

% Wavefunction ansatz
spin=0;
psi=@(a1,a2) real(a1*exp(1i*spin*tt-(rr/a2).^2).*(rr.^spin));
% Cost function
f=@(a1,a2) lag(psi(a1,a2));

% Initial guess
a=[(0.5767/2)^((spin+1)/2); 3.456];
a=[sqrt(2); 2];
vars=num2cell(a);

setlatex();
figure(1);
hp=surf(xx,yy,reshape(abs(psi(vars{:})),[n,n]));
shading interp;
axis tight;
colormap(magma(256));
view(2);
colorbar();

% Newton-Raphson for gradient
g=agrad(f,length(a));
J=ahess(f,length(a));
iter=0;
y=ones(size(a));
tol=1e-12;
while( norm(y)>tol && iter<60 )
    vars=num2cell(a);
    y=g(vars{:});
    a=a-J(vars{:})\y;
    iter=iter+1;
    
    set(hp,'ZData',abs(psi(vars{:})));
    E=real(f(vars{:}));
    title(num2str(E,'$E = %f$'))
    drawnow;
end

display(iter);
display(y);
display(a);
end