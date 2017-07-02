function [] = heatExp( N )
% Linear PDE with mixed boundary conditions
%
% u_t = u_xx
% u + c*du/dn = constant
%
% Method of lines: discrete x, continous t
% u_t = D2*u -> u(t) = expm(t*D2)*u(0)
%
% Constrained evolution is achieved via projection method

c=1;
dt=0.001;

[D,x]=chebD(N);
A=D*D;

a=[1,-1];
b=[c, c];

I=eye(size(A,1));
B=diag(a)*I([1,end],:)+diag(b)*D([1,end],:); 
P=I-B'/(B*B')*B;
R=null(A);
T=R/(R'*R)*R';
Q=T+expm(dt*P*A)*(I-T);

u=exp(-10*x.^2)-exp(-10);

figure(1);
h=plot(x, u);
xlim([-1,1]);
ylim([-4,4]);

nframes=1000;
for i=1:nframes
    u=Q*u;
    set(h, 'YData', u);
    %disp(B*u);
    drawnow;
end
end