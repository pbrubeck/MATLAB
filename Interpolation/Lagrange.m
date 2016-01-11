function [p] = Lagrange(x, y, t)
% Evaluates the Lagrange interpolation polynomial.
n=length(x); 
m=length(t);
X=repmat(x(:).',[n,1]);
w=prod(X-X.'+eye(n));
r=y./w;
D=repmat(t(:).',[n,1])-repmat(x(:),[1,m]);
p=prod(D).*(r*(1./D));
end