function [] = brachistochrone(n)
global D x w y1 y2 h;
L=20;

[~,w]=ClenshawCurtis(-1,1,n); w=w(:);
[D,x]=chebD(n);

y1=0;
y2=2;
y0=y1+(y2-y1)*(x(2:end-1)+1)/2;
h=max(y1,y2);
x=L/2*(x+1); D=2/L*D;

opts = optimoptions(@fmincon,'Algorithm','interior-point',...
    'GradObj','off','MaxFunEvals',10000,'PlotFcns',@simpleplot);
lb=(min(y1,y2)-L)*ones(n-2,1);
ub=(max(y1,y2)+L)*ones(n-2,1);

fmincon(@gradfnc,y0,[],[],[],[],lb,ub,[],opts);
end

function [f]=gradfnc(y)
global D w y1 y2 h;
y=h-[y2;y(:);y1];
dy=D*y;
v0=1;

f=w'*(sqrt(1+dy.^2)./(v0+sqrt(y)));
end

function [stop]=simpleplot(y,optimValues,state,varargin)
global x y1 y2;
y=[y2;y(:);y1];
plot(x,y); hold on;
stem(x,y); hold off;
title(num2str(optimValues.fval,'Time = %f'));
stop = false;
end
