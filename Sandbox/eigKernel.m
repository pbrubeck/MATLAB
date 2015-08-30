n=256; p=0;
[x, w]=GaussLegendre(-4, 4, n);


x2 = x.^2;
b=exp(0.5i*x2);
K = (besselj(p, x*x').*(b*b.'))*diag(x.*exp(-x2))*diag(w);
[V,D]=eig(K,'nobalance');
lambda=diag(D);
norm(K-V*D/V, 'fro')


m=6;
disp(lambda(1:m));
figure(1);
plot(x,abs(V(:,1:m)));
figure(2);
plot(x,angle(V(:,1:m)));

mesh(x,x,real(K))