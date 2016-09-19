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
B=diag(x);
T=-1/2*(diag(x)*D*D+2*D);
V=-eye(N);
H=T+V;
[S,E]=eig(H(1:end-1,1:end-1), B(1:end-1,1:end-1)); 
[E,id]=sort(diag(E));
E=E0*E(1:k);

% Normalize eigenfunction
psi=zeros(N,k);
psi(1:end-1,:)=S(:,id(1:k));
psi=normc(psi, w.*r.^2);
psi=bsxfun(@times, psi, sign(psi(1,:)));

% Calculate Bohr radius
u=r.*psi(:,1);
uu=abs(u).^2;
[Pmax, imax]=max(uu);
id=imax-2:imax+3;
p=polyfit(r(id), uu(id),5);
dp=polyder(p);
rs=roots(dp);
rs=rs(abs(imag(rs))<eps & real(rs)>0);
rbohr=rs(find(min(abs(r(imax)-rs))));


% Plot
plot(r, psi);
xlim([0,30]);
legend([num2str((0:k-1)','E_{%d}'), num2str(E,'=%f')]);
end