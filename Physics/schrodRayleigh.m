function [] = schrodRayleigh(n, m)

x=linspace(0,1,n+2); x=x(:);
V=2000*(x-0.5).^2;
%V=2000./abs(x-0.5);
%V=(x<0.7).*(-500*x);
%V=400*((x<0.6)-(x<0.4));
M=diag((pi*(1:n)).^2)+2/(n+1)*dst(dst(diag(V(2:end-1)))')';
[S, E]=eigs(M, m, 'sm'); [E, idx]=sort(diag(E)); S=S(:,idx);
S=sqrt(n+1)*dst(S);
S=S/diag(sqrt(diag(S'*S)));
Psi=zeros(n+2, m);
Psi(2:end-1,:)=S;

figure(1); plot(x, Psi);
v=[1:m; E']; title(sprintf('E%i=%f ', v(:)));

a=exp(2i*pi*rand(m,1))./(E.^2); a=a./sqrt(a'*a);
figure(2);
h=plot(x, zeros(size(x)));
ylim([0, 0.1]);
for t=0:0.0001:10
    u=Psi*(a.*exp(-1i*E*t));
    h.YData=abs(u).^2;
    drawnow;
end
end