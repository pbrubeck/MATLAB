function [] = jumpStiff(n,xi)

k=1;
L=1;

N = (n-1)*k+1; % Total number of gridpoints
mask=(mod(0:(n-1)*(k-1),n-1)==0)';

% Element Diff matrix, nodes and quadrature
[D0,x0,w0]=legD(n);
V0=VandermondeLeg(x0);
M0=inv(V0*V0'); % Element mass
K0=diag(w0)*D0; % Element stiffness

% Element boundaries
xe=linspace(-L,L,k+1)';
J0=diff(xe)/2; % Jacobian
x1=kron(J0, x0)+kron((xe(1:end-1)+xe(2:end))/2, ones(n,1));
x1(n+1:n:end)=[];

% Spectral element assembly
J1=kron(J0,ones(n-1,1));
J1=J1(1:numel(mask));
M1=conv2(diag(mask.*J1), M0);
K1=conv2(diag(mask    ), K0);


MJ=jumpForce(xi,xe,x0,x1,M0,-eye(n),zeros(1,n),-1);
KJ=jumpForce(xi,xe,x0,x1,K0,-eye(n),zeros(1,n), 0);

display(MJ);
display(KJ);
end