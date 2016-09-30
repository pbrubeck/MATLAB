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

lb=(min(y1,y2)-L)*ones(n-2,1);
ub=(max(y1,y2)+L)*ones(n-2,1);
path=y1+(y2-y1)*(x(2:end-1)+1)/2;
[path,t]=ga(@timeFun,n-2,[],[],[],[],lb,ub,[],[],opts);
disp(t);
[path,t]=fmincon(@timeFun,path,[],[],[],[],lb,ub);
disp(t);

y=[y2;path(:);y1];
plot(x,y);
end

function t=timeFun(y)
global D w y1 y2 h;
y=[y2;y(:);y1];
t=w*sqrt((1+(D*y).^2)./(h-y));
if ~isreal(t)
    t=inf;
end
end
