% Sizes
n=32;    % Chebyshev grid
nj=16;    % Derivative jumps enforced
Nq=1024; % Interpolation grid

% Domain
a=0;     % Lower limit
b=5;     % Upper limit
xi=10/3; % Discontinuity
[D,x0]=chebD(n);
x=(b-a)/2*(x0+1)+a;
xx=linspace(a,b,Nq)';

% Piecewise Function
coef=[0,0,1];
f1=@(x) real(LegendreQ(coef,xi))*LegendreP(coef,x);
f2=@(x) LegendreP(coef,xi)*real(LegendreQ(coef,x));
fun=@(x) (x<xi).*f1(x)+(x>=xi).*f2(x);

% Get jumps
z0=ainit(xi,nj-1);
y1=f1(z0);
y2=f2(z0);
jumps=zeros(nj,1);
for r=1:nj
    jumps(r)=y2{r-1}-y1{r-1};
end

% Interpolate
[s1,s2]=piecewiseLagrange(x,xi,jumps);
P=interpcheb(eye(n),linspace(-1,1,Nq)');
y=fun(x);
yy=zeros(size(xx));
yy(xx<=xi)=P(xx<=xi,:)*(y+s1);
yy(xx>=xi)=P(xx>=xi,:)*(y+s2);

figure(1);
plot(x,y,'.k', xx,yy,'r');
title('Interpolation');

figure(2);
plot(xx,fun(xx)-yy);
title('Error');