function [y] = LegendreP(varargin)
% LegendreP(a,x) Evaluates the Legendre series given by the coeficients a(n)
% LegendreP(a,m,x) Associated Legendre Functions
a=varargin{1};
x=varargin{end};
if(nargin==2)
    m=0;
elseif(nargin==3)
    m=varargin{2};
end
n=length(a);
y=(n>m+1)*a(n);
yy=zeros(size(x));
for k=n-2:-1:m+1
    temp=y;
    y=a(k+1)+(2*k+1)/(k-m+1)*x.*y-(k+m+1)/(k-m+2)*yy;
    yy=temp;
end
Pmm=(-1)^m*prod(1:2:2*m-1)*(1-x.^2).^(m/2);
y=Pmm.*(a(m+1)+(2*m+1)*(x.*y-yy/2));
end