function [] = HermiteODE(n, m)
[D, x, w]=hermD(n+1);
H=-D*D+diag(x.^2);
[S, lambda]=eigs(H(2:end-1,2:end-1), m, 'sm');
lambda=diag(lambda);
Psi=zeros(n+2, m);
Psi(2:end-1,:)=S;

%Normalization
w=w.*exp(x'.^2);
Psi=bsxfun(@rdivide, Psi, sqrt(w*(Psi.^2)));

plot(x, Psi(:,end)); 
disp(lambda);
end