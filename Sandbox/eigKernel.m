n=512; p=0; a=-4; b=4;
[x,w]=gauleg(a,b,n);

x2 = x.^2;
r=exp(0.5i*x2);
K=diag(r)*besselj(p, x'*x)*diag(w.*r.*x.*exp(-x2));

[V,D]=eig(K);
lambda=diag(D);

m=3;
figure(1);
disp(lambda(1:m));
plot(x,abs(V(:,1:m)));
figure(2);
plot(x,angle(V(:,1:m)));