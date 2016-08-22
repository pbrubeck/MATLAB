function [E, rbohr] = schrodHydrogen(N, k)
% Schrodinger Equation for the Hydrogen atom
E0=54.4228;  % Electronvolts
r0=0.264589; % Angstroms

% Differentiation matrix
[D,x,w]=lagD(N);
r=r0*x;
w=r0*w(:);

% Solve eigensystem
T=-D*diag(x.^2)*D;
V=-diag(x);
H=T+V;
B=diag(x.^2);
[S,E]=eig(H(1:end-1,1:end-1), B(1:end-1,1:end-1));
[E,id]=sort(diag(E));
E=E0*E(1:k);

% Normalize eigenfunction
Psi=zeros(N,k);
Psi(1:end-1,:)=S(:,id(1:k));
Psi=Psi/(diag(sqrt(diag(Psi'*diag(w)*Psi))));
Psi=bsxfun(@times, Psi, sign(Psi(1,:)));

% Calculate Bohr radius
rbohr=Psi(:,1)'*diag(2*r.*w)*Psi(:,1);

% Plot
plot(r, Psi);
xlim([0,20]);
legend([num2str((1:k)','E_{%d}'), num2str(E,'=%f')]);
end

