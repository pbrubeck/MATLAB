% File GreatLakes.m
% Honors ODE final assesment

% Data
V=[2900,0,0,0,0;0,1180,0,0,0;0,0,850,0,0;0,0,0,116,0;0,0,0,0,393];
G=[-15,0,0,0,0;0,-38,0,0,0;15,38,-68,0,0;0,0,68,-85,0;0,0,0,85,-99];
x0=V*[1; 1; 1; 1; 1];

A=G/V;
[P,D]=eig(A);
d=diag(D);

PC=P*diag(P\x0);
PCD=PC*D;

t=linspace(0, 100, 100);
xt=PC*exp(d*t);
rt=diag(x0)\xt;

rect=diag(V)*ones(size(t))+diag(G)*t;
approx=V*exp(V\diag(G)*t);
figure(1)
plot(t,rect,'LineWidth',2)
legend('Superior','Michigan','Huron','Erie','Ontario');
axis([0,100,0,3000]);
figure(2)
plot(t,approx,'LineWidth',2)
legend('Superior','Michigan','Huron','Erie','Ontario');
figure(3)
hold on
plot(t,ones(1,5)*approx,'LineWidth',2)
plot([0 100],[5167.05 5167.05],'LineWidth',2)
plot([0 100],[2719.5 2719.5],'LineWidth',2)
hold off

figure(4);
plot(t,xt,'LineWidth',2);
legend('Superior','Michigan','Huron','Erie','Ontario');
figure(5);
plot(t,rt,'LineWidth',2);
legend('Superior','Michigan','Huron','Erie','Ontario');


k=[0.95, 0.5];
T=zeros([length(x0), length(k)]);
for i=1:length(x0)
    B=[PC(i,:); PCD(i,:)];
    for j=1:length(k)
        tx=0;
        it=0;
        err=1;
        xk=x0(i)*k(j);
        while(abs(err)>1E-15 && it<20)
            f=B*exp(d*tx);
            err=(f(1)-xk)/f(2);
            tx=tx-err;
            it=it+1;
        end
        T(i,j)=tx;
    end
end
disp(k);
disp(T);