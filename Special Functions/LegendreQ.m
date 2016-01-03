function [y] = LegendreQ(varargin)
% LegendreQ(a,x) Evaluates the Legendre series given by the coeficients a(n)
% LegendreQ(a,m,x) Associated Legendre Polynomials
a=varargin{1};
x=varargin{end};
if(nargin==2)
    m=0;
elseif(nargin==3)
    m=varargin{2};
end
n=length(a);
y=(n>m+1)*a(n); yy=zeros(size(x));
for k=n-2:-1:m+1
    temp=y;
    y=a(k+1)+(2*k+1)/(k-m+1)*x.*y-(k+m+1)/(k-m+2)*yy;
    yy=temp;
end
tau=atanh(x);
y=tau*a(m+1)+(x.*tau-1).*y-tau.*yy/2;
end