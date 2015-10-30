function [] = Schrodinger(n, m)
[D, x]=chebD(n+2); x=(1+x)/2;
D2=4*D^2; D2=D2(2:end-1,2:end-1);
V=2000*(x-0.5).^2;
%V=(x<0.7).*(-500*x);
H=-D2+diag(V(2:end-1));

[S, E]=eig(H, 'vector');
[foo, idx]=esort(-E); E=E(idx);

y=linspace(0,1,n+2); y=y(:);
Psi1=zeros(n+2, m);
Psi1(2:end-1,:)=S(:,idx(1:m));
Psi1=normc(interp1(x, Psi1, y, 'spline'));
figure(1); plot(y, Psi1);
v=[1:m; E(1:m)'];
title(sprintf('E%i=%f ', v(:)));



x=linspace(0,1,n+2); x=x(:);
V=2000*(x-0.5).^2;
%V=(x<0.7).*(-500*x);
M=diag((pi*(1:n)).^2)+2/(n+1)*dst(dst(diag(V(2:end-1)))')';

[S, E]=eig(M, 'vector');
[foo, idx]=esort(-E); E=E(idx);

Psi2=zeros(n+2, m);
Psi2(2:end-1,:)=normc(dst(S(:,idx(1:m))));
figure(2); plot(x, Psi2);
v=[1:m; E(1:m)'];
title(sprintf('E%i=%f ', v(:)));

disp(abs(diag(Psi1'*Psi2)))
end
