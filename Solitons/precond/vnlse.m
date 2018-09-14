function [E] = vnlse(n)
% Stiffness matrix
[D,x,w]=legD(n);
hx=4;
x=hx*(x+1);
D=D/hx;
w=x.*w*hx;
K=D'*diag(w)*D;
M=diag(w);
lam1=1/2;
lam2=1;
H=lam1*M+lam2*K;

% Nonlinear term
nl=@(psi2) -psi2.^2/2;
% Lagrangian
lag=@(psi) (psi'*H*psi + w'*nl(psi.*conj(psi)))/2;

% Wavefunction ansatz
psi=@(a1,a2) a1*exp(-(x/a2).^2);
% Cost function
f=@(a1,a2) lag(psi(a1,a2));

% Initial guess
a=[1.5;1.5];
vars=num2cell(a);

setlatex();
figure(1);
hp=plot(x,abs(psi(vars{:})),'-ok');

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
    
    set(hp,'YData',abs(psi(vars{:})));
    E=f(vars{:});
    title(num2str(E,'$E = %f$'))
    drawnow;
end

display(iter);
display(y);
display(a);
end