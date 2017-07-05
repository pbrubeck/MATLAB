function [] = burgersDG(k,p)
xe=linspace(-1,1,k+1)'; % Element boundaries
[Dx,xn,wn]=legD(p); % Diff matrix, nodes and quadrature
J=diff(xe)/2; % Jacobian
x=kron(J, xn)+kron((xe(1:end-1)+xe(2:end))/2,ones(p,1));

% Mass matrix
V=VandermondeLeg(xn);
Minv=V*V';

% Stiffness matrix
K=Dx'*diag(wn);

mass=@(u) reshape(Minv*reshape(u,p,k)./J',[],1);
stiff=@(u) reshape(K*reshape(u,p,k),[],1);

% Lax-Friedrichs numerical flux
function [F]=LaxFriedrichs(u)
    uu=reshape(u,p,k);
    u0  =(uu(end,[end,1:end-1])+uu(1,1:end))/2;
    jump=(uu(end,[end,1:end-1])-uu(1,1:end))/2;
    f=u0.^2/2+max(abs(uu)).*jump;
    F=zeros(p,k);
    F(end,[end,1:end-1])=f;
    F(1,:)=-f;
    F=reshape(F,[],1);
end

function du=partialT(u)
    F=LaxFriedrichs(u);
    du=mass(stiff(u.^2/2)-F);
end

function u=solveRK4(u,dt)
    k1=dt*partialT(u);
    k2=dt*partialT(u+k1/2);
    k3=dt*partialT(u+k2/2);
    k4=dt*partialT(u+k3);
    u=u+(k1+2*k2+2*k3+k4)/6;
end

% Initial condition
u=(sin(pi*x).^2).*(x<0);

figure(1);
h1=plot(x,u);
axis manual;

% Time propagation
dt=0.0005;
T=10;
nframes=ceil(T/dt);
for i=1:nframes
    u=solveRK4(u,dt);
    set(h1,'YData',u);
    drawnow;
end
end