function [] = eulerDG(k,p)

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
    
end

function du=partialT(u)
    [F,G]=eulerflux(u);
    f=LaxFriedrichs(u);
    du=mass(stiff(f)-f);
end

function u=solveRK4(u,dt)
    k1=dt*partialT(u);
    k2=dt*partialT(u+k1/2);
    k3=dt*partialT(u+k2/2);
    k4=dt*partialT(u+k3);
    u=u+(k1+2*k2+2*k3+k4)/6;
end


end


function [F,G]=eulerflux(q)
    rho=q(:,:,1);
    rhou=q(:,:,2);
    rhov=q(:,:,3);
    E=q(:,:,4);
    u=rhou./rho;
    v=rhov./rho;
    
    p=(2/5)*(E-(rhou.^2+rhov.^2)./(2*rho));

    F(:,:,1)=rhou;
    F(:,:,2)=rhou.*u+p;
    F(:,:,3)=rhou.*v;
    F(:,:,4)=u.*(E+p);
    
    G(:,:,1)=rhov;
    G(:,:,2)=rhov.*u;
    G(:,:,3)=rhov.*v+p;
    G(:,:,4)=v.*(E+p);
end
