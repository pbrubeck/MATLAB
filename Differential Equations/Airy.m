function [u] = Airy(m, n)
% Solves Airy's ordinary differential equation on the interval containing
% the first m zeros
[D, x]=chebD(n);
D2=D^2; D2=D2(2:n-1, 2:n-1);
[V,lambda]=eig(D2, diag(x(2:n-1)), 'vector');
ii=find(lambda>0); 
lambda=lambda(ii); V=V(:,ii);
[~, ii]=sort(lambda); ii=ii(m);
lambda=lambda(ii);
u=[0; V(:, ii); 0];
mid=(u(ceil(end/2))+u(floor(end/2+1)))/2;
u=(airy(0)/mid)*u;
plot(x, u); title(sprintf('lambda = %f', lambda));
end