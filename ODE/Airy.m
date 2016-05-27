function [u] = Airy(m, n)
% Solves Airy's ordinary differential equation on the interval containing
% the first m zeros
[D, x]=chebD(n);
D2=D^2; D2=D2(2:n-1, 2:n-1);
[V,lambda]=eigs(D2, diag(x(2:end-1)), m, 'sm');
lambda=diag(lambda);
[~, ii]=sort(lambda); ii=ii(end); 
lambda=lambda(ii);
u=[0; V(:, ii); 0];
mid=(u(ceil(end/2))+u(floor(end/2+1)))/2;
u=(airy(0)/mid)*u;
figure(1); plot(x, u);
title(sprintf('lambda = %f', lambda));
end