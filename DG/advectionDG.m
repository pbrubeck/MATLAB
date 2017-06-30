function [] = advectionDG(k,p)
xe=linspace(-1,1,k+1)'; % Element boundaries
[Dx,xn,wn]=legD(p); % Diff matrix, nodes and quadratures
J=diff(xe)/2;
x=kron(J, xn)+kron((xe(1:end-1)+xe(2:end))/2,ones(p,1));

% Full operators
W=kron(J, wn);
D=kron(diag(1./J), Dx);

% Matrix-free operators
w=@(x) wn'*reshape(x,p,k)*J;
d=@(x) reshape(Dx*reshape(x,p,k)./J', size(x));

% Advection Numerical Flux
c=-1; % velocity
s=1; % s=0 average, s=1 upwind
F=zeros(p,p+2);
F0=[1,1;-1,-1]/2;
F0=c*F0+abs(c)*s*F0';
F(1,1:2)=F0(1,:);
F(end,end:-1:end-1)=F0(end,:);
disp(F)

% Mass matrix
%M=diag(wn)*J(1);
M=J(1)*[2/3 1/3; 1/3 2/3];

% Stiffness matrix
K=zeros(p,p+2);
%K(:,2:end-1)=Dx'*diag(wn);
K(:,2:end-1)=[-1/2 -1/2; 1/2 1/2];

% Stencil
S=M\(K+F);

A = zeros(k*p);
for i=1+size(S,1) : size(S,1) : size(A,2)-size(S,1)
    A(i:i-1+size(S,1),i-1:i-2+size(S,2))=S;
end
A(1:size(S,1),[end, 1:size(S,2)-1])=S;
A(end+1-size(S,1):end,[end+2-size(S,2):end, 1])=S;

u=sin(2*pi*x);
h1=plot(x,u);

CFL=0.1;
dt=CFL*max(J)/abs(c);

% Time propagator
Q=eye(k*p)+dt*A;

T=3;
nframes=ceil(T/dt);
for i=1:nframes
    u=Q*u;
    set(h1,'YData',u);
    drawnow;
end

end

