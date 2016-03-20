function [] = LegendreODE(n, m)
[D, x]=chebD(n+2);
P=D*(diag(1-x.^2)*D);
P=P(2:end-1,2:end-1);
[S, E]=eigs(P, m, 'sm'); 
[E, idx]=sort(diag(E));
S=S(:,idx);
Psi=zeros(n+2, m);
Psi(2:end-1,:)=S;
plot(x, Psi); disp(E);
end