function [E, rbohr] = schrodHydrogen(N, k)
% Schrodinger Equation for the Hydrogen atom
% Hatree units
E0=27.21138505;   % Electronvolts
r0=0.52917721067; % Angstroms

% Differentiation matrix
[D,x,w]=lagD(N);
r=r0*x;
w=r0*w;

% Solve eigensystem
T=-1/2*D*diag(x.^2)*D;
V=-diag(x);
H=T+V;
B=diag(x.^2);
[S,E]=eig(H(1:end-1,1:end-1), B(1:end-1,1:end-1)); 
[E,id]=sort(diag(E));
E=E0*E(1:k);

% Normalize eigenfunction
Psi=zeros(N,k);
Psi(1:end-1,:)=S(:,id(1:k));
Psi=normc(Psi, w.*r.^2);
Psi=bsxfun(@times, Psi, sign(Psi(1,:)));

% Calculate Bohr radius
u=r.*Psi(:,1);
[Pmax, imax]=max(abs(u).^2);
rq=linspace(r(imax-1), r(imax+1), 1024);
uq=interp1(r,u,rq,'spline');
[Pmax, imax]=max(abs(uq).^2);
rbohr=rq(imax);

% Plot
plot(r, Psi);
xlim([0,20]);
legend([num2str((0:k-1)','E_{%d}'), num2str(E,'=%f')]);
end