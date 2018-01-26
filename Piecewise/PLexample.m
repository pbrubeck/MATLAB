% Sizes
n=32;    % Chebyshev grid
nj=16;    % Derivative jumps enforced
Nq=1024; % Interpolation grid

% Domain
a=0;     % Lower limit
b=5;     % Upper limit
xi=10/3; % Discontinuity
[D,x]=chebD(n);
x=(b-a)/2*(x+1)+a;
xx=linspace(a,b,Nq)';

% Piecewise Function
coef=[0,0,1];
f1=@(x) real(LegendreQ(coef,xi))*LegendreP(coef,x);
f2=@(x) LegendreP(coef,xi)*real(LegendreQ(coef,x));
fun=@(x) (x<xi).*f1(x)+(x>=xi).*f2(x);

% Get jumps
x0=ainit(xi,nj-1);
y1=f1(x0);
y2=f2(x0);
jumps=zeros(nj,1);
for r=1:nj
    jumps(r)=y2{r-1}-y1{r-1};
end

% Interpolate
s=piecewiseLagrange(x,jumps);
P=interpcheb(eye(n),linspace(-1,1,Nq));

figure(3)
imagesc((xx'>=xi).*(xi>x)-(xx'<xi).*(xi<x))

u=fun(x);
y=P*u+diag(P*s(xx',xi));

figure(1);
plot(xx,y, xx,fun(xx));
title('Interpolation');

figure(2);
plot(xx,fun(xx)-y);
title('Error');