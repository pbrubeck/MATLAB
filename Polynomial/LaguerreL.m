function [y] = LaguerreL(varargin)
% LaguerreL(c,x) Evaluates the Laguerre series given by the coefficients c(n)
% LaguerreL(c,a,x) Associated Laguerre Polynomials
c=varargin{1};
x=varargin{end};
if(nargin==2)
    a=0;
elseif(nargin==3)
    a=varargin{2};
end
n=length(c);
y=(n>1)*c(n); yy=zeros(size(x));
for k=n-2:-1:1
    temp=y;
    y=c(k+1)+(2*k+1+a-x)/(k+1).*y-(k+1+a)/(k+2)*yy;
    yy=temp;
end
y=c(1)+(1+a-x).*y-(1+a)/2*yy;
end