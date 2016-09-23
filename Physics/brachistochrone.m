function [] = brachistochrone(n)
global D w y1 y2 h;
L=12;

[~,w]=ClenshawCurtis(-1,1,n-1);
[D,x]=chebD(n);
x=L/2*(x+1); D=2/L*D;

y1=2;
y2=0;
h=max(y1,y2)+L;
opts=[];
opts.FunctionTolerance=1e-13;

[path,t]=ga(@timeFun,n-2,[],[],[],[],min(y1,y2)-L,max(y1,y2)+L,[],[],opts);
y=[y2;path(:);y1];
plot(x,y);
disp(t);
end

function t=timeFun(y)
global D w y1 y2 h;
y=[y2;y(:);y1];
t=w*sqrt((1+(D*y).^2)./(h-y));
if ~isreal(t)
    t=inf;
end
end
