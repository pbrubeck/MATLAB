function []=pendulum(th0, w0, F, n)
h=60/(n-1);
u=[th0, w0];
ut=zeros(2, n); ut(:,1)=u;
ti=0;
for i=1:n-1
    k1=h*timeD(ti, u, F);
    k2=h*timeD(ti+h/2, u+k1/2, F);
    k3=h*timeD(ti+h/2, u+k2/2, F);
    k4=h*timeD(ti+h, u+k3, F);
    u=u+(k1+2*k2+2*k3+k4)/6;
    ti=ti+h;
    ut(:, i+1)=u; 
end
t=linspace(0,60,n);
figure(1);
plot(t, ut);
end

function v=timeD(t, u, F)
v=[u(2), -sin(u(1))-0.5*u(2)+F*sin(2*t/3)];
end

