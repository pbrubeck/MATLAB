function [] = HermiteODE(n, m)
L=8; 
L2=L*L;
[D, x]=chebD(n+2); x=L*x;
D2=D^2/L2; 
D2=D2(2:end-1,2:end-1);
V=x.^2;
H=-D2+diag(V(2:end-1));
[S, E]=eigs(H, m, 'sm');
[E, idx]=sort(diag(E));
S=S(:,idx);
S=bsxfun(@times, S, sign(S(1,:)));
Psi=zeros(n+2, m);
Psi(2:end-1,:)=S;

% Normalization
w=(pi*L/(n+1))*sqrt(1-(x/L).^2);
Psi=bsxfun(@rdivide, Psi, sqrt(w'*(Psi.^2)));

plot(x, Psi); 
disp(E);
end