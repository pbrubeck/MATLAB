function [] = nonlinearODE(n)
[D,x]=chebD(n);
D2=D*D;
u1=1;
u2=2;

bc=D2(2:end-1,[1 end])*[u2;u1];
u=u1+(u2-u1)*(x+1)/2;

i=0; err=1; tol=1e-10;
while err>tol
    unew=D2(2:end-1,2:end-1)\(exp(u(2:end-1))-bc);
    err=norm(unew-u(2:end-1),inf);
    u(2:end-1)=unew;
    i=i+1;
end
fprintf('%d iterations\n',i);
plot(x,u);
end

