% File GreatLakes.m
% Honors ODE final assesment

% Data
V=[2900,0,0,0,0;0,1180,0,0,0;0,0,850,0,0;0,0,0,116,0;0,0,0,0,393];
C=[-15,0,0,0,0;0,-38,0,0,0;15,38,-68,0,0;0,0,68,-85,0;0,0,0,85,-99];

A=C/V;
[P,D]=eig(A);
d=diag(D);

x0=V*[1 1 1 1 1]';
b=P\x0;

t=linspace(0, 100, 1000);
xt=zeros(size(A,1), size(t,2));
for i=1:size(t,2)
    xt(:,i)=P*diag(exp(d*t(i)))*b;  
end



rt=diag(x0)\xt;
disp(d);

figure(1);
plot(t,xt);
figure(2);
plot(t,rt);