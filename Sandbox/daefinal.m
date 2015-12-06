y=[2.83 2.1 1.71
2.63 2.31 1.55
1.84 2.50 1.70
2.70 2.98 1.76
2.15 2.65 1.89
2.91 2.81 2.10
2.46 3.16 1.38
2.31 2.36 1.60
2.33 2.63 1.28];

y_hat=ones(9,1)*sum(A)/9;
err=(y-y_hat)/(var(y(:)));
h1=figure(1);
scatter(y_hat(:),err(:),'o','k','filled');
xlabel('Estimador $\hat{y}$','interpreter','latex');
ylabel('Error $\epsilon$','interpreter','latex');

z=norminv(((1:27)'-0.5)/27,0,1);
err=sort(err(:));

X = [ones(length(err),1) err];
b = X\z;
zCalc=X*b;
R2=1-sum((z-zCalc).^2)/sum((z-mean(z)).^2);
h2=figure(2); hold on;
scatter(err,z,'o','k','filled');
plot(err, zCalc, 'r'); 
text(-1, 1, sprintf('$R^2$ = %f', R2), 'interpreter','latex');
hold off;
xlabel('Error $\epsilon$','interpreter','latex');
ylabel('Percentil $z$','interpreter','latex');

