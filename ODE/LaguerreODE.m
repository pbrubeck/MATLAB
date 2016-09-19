function [E] = LaguerreODE( N, k, a )
[D,x,w]=lagD(N);

A=-(eye(N)+D)*diag(x.^(a+1))*D;
B=diag(x.^a);

[V,E]=eig(A(1:end-1,1:end-1),B(1:end-1,1:end-1),'vector');
[E,id]=sort(E);

psi=zeros(N,k);
psi(1:end-1,:)=V(:,id(1:k));
psi=bsxfun(@times, conj(psi(1,:)), psi);
psi=normc(psi, w.*(x.^a));
psi=bsxfun(@times, sqrt(gamma(1+a+(0:k-1))/factorial(0:k-1)), psi);
plot(x,real(psi)); xlim([0, 10]);
end

