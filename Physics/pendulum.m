function [ ] = pendulum( th0, w0, F, n )
h=60/(n-1);
g=9.8;
q=0.5;
L=9.8;
W=2/3;

u=[th0, w0];
ut=zeros(2, n); ut(:,1)=u;
ti=0;
for i=1:n-1
    k1=h*timeD(ti, u, g, L, q, F, W);
    k2=h*timeD(ti+h/2, u+k1/2, g, L, q, F, W);
    k3=h*timeD(ti+h/2, u+k2/2, g, L, q, F, W);
    k4=h*timeD(ti+h, u+k3, g, L, q, F, W);
    u=u+(k1+2*k2+2*k3+k4)/6;
    ti=ti+h;
    ut(:, i+1)=u; 
end
t=linspace(0,60,n);
E=g*L*(1-cos(ut(1,:)))+(L*ut(2,:)).^2/2;
figure(1);
plot(t, ut);
figure(2);
plot(t, E);
end

function v=timeD(t, u, g, L, q, F, W)
v=[u(2), -g/L*sin(u(1))-q*u(2)+F*sin(W*t)];
end

