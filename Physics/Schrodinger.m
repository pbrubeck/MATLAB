function [] = Schrodinger(n, m)
tic();
[D, x]=chebD(n+2); x=(1+x)/2;
D2=4*D^2; D2=D2(2:end-1,2:end-1);
%V=2000./abs(x-0.5);
%V=2000*(x-0.5).^2;
%V=(x<0.7).*(-500*x);
V=400*((x<0.6)-(x<0.4));
H=-D2+diag(V(2:end-1));
[S, E]=eigs(H, m, 'sm'); [E, idx]=sort(diag(E)); S=S(:,idx);
Psi1=zeros(n+2, m);
Psi1(2:end-1,:)=S;
toc()

y=linspace(0,1,n+2); y=y(:);
Psi1=sqrt(n+1)*normc(interp1(x, Psi1, y, 'spline'));
figure(1); plot(y, Psi1);
v=[1:m; E']; title(sprintf('E%i=%f ', v(:)));


tic();
x=linspace(0,1,n+2); x=x(:);
%V=2000./abs(x-0.5);
%V=2000*(x-0.5).^2;
%V=(x<0.7).*(-500*x);
V=400*((x<0.6)-(x<0.4));
M=diag((pi*(1:n)).^2)+2/(n+1)*dst(dst(diag(V(2:end-1)))')';
[S, E]=eigs(M, m, 'sm'); [E, idx]=sort(diag(E)); S=S(:,idx);
Psi2=zeros(n+2, m);
Psi2(2:end-1,:)=sqrt(n+1)*normc(dst(S));
toc()

figure(2); plot(x, Psi2);
v=[1:m; E']; title(sprintf('E%i=%f ', v(:)));

a=normc(1./E);
figure(3);
h=plot(x, zeros(size(x)));
ylim([0, 5]);
for t=0:0.0005:10
    h.YData=abs(Psi2*(a.*exp(-1i*E*t))).^2;
    drawnow();
end
end