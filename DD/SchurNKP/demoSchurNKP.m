function [lam] = demoSchurNKP(m, shape, ref)
% m = Number of gridpoints in 1D per quadrilateral element
% shape = 1, Equilateral triangle
% shape = 2, Right triangle
% shape = 3, Isoceles triangle
% shape = 4, Scalene triangle
% shape = 5, L shaped membrane
% shape = 6, Unit disc
% shape = 7, Annulus
% shape = 8, Newtonian binary stars
% ref = Level of refinement (0 = no refinement)

function F=myrhs(x,y)
    F=1+0*x+0*y;
end
L=3;
R1=1;  R2=1;
z1=L;  z2=-L;
m1=1;  m2=1;
function F=binrhs(z,rho)
    r1=hypot(rho,z-z1);
    r2=hypot(rho,z-z2);
    F1=-m1*(r1<=R1).*sinc(pi*r1/R1);
    F2=-m2*(r2<=R2).*sinc(pi*r2/R2);
    F=F1+F2;
end

RHS=@myrhs;
metric = @(x,y) 1+0*x+0*y;

% Triangle Vertices
tri=zeros(3,4);
tri(:,1)=[1i;exp(1i*pi*7/6);exp(-1i*pi*1/6)]; % Equilateral
tri(:,2)=2/(2-sqrt(2))*[1i;0;1];              % Right angle
tri(:,3)=[2i;-1;1];                           % Isoceles
tri(:,4)=[-2+3i;0;2];                         % Scalene

if(shape<=size(tri,2))
    L=abs(tri([3,1,2],shape)-tri([2,3,1],shape)); % Sides
    V=eye(3)+diag((sum(L)/2-L)./L([2,3,1]))*[-1,0,1; 1,-1,0; 0,1,-1];
    z0=zeros(7,1);
    z0([1,2,3])=tri(:,shape);       % Vertices
    z0([4,5,6])=V*tri(:,shape);     % Touch points
    z0(7)=(L'*tri(:,shape))/sum(L); % Incenter
    quads=[7,5,4,1; 7,6,5,2; 7,4,6,3];
    curv=zeros(size(quads)); curv(:)=inf;
else
    switch(shape-size(tri,2))
        case 1
            [z0,quads,curv]=Lgrid();
        case 2
            [z0,quads,curv]=coregrid(1);
            metric = @(z,rho) abs(rho);
        case 3
            [z0,quads,curv]=annulusgrid(1,2);
            metric = @(z,rho) abs(rho);
        case 4
            [z0,quads,curv]=bingrid(R1,R2,L);
            metric = @(z,rho) abs(rho);
            RHS=@binrhs;
    end
end

% Mesh refinement
for j=1:ref
    [~, adj, bnd] = meshtopo(quads);
    [z0,quads,curv]=quadmeshrefine(z0,quads,curv,adj,bnd);
end

k=32;
[lam, uu, X, Y, relres, pcalls] = meshSchurNKP(m, z0, quads, curv, k, RHS, metric);
delete(findall(gcf,'Type','light'));
% shading faceted;
display(relres);
display(pcalls);

if shape==6
    r=hypot(X,Y);
    uex=(1-r.^2)/6;
    
    err=norm(uu(:)-uex(:),'inf')/norm(uex(:),'inf');
    disp(err);
    
    figure(2);
    for j=1:size(quads,1)
        surf(X(:,:,j), Y(:,:,j), uu(:,:,j)-uex(:,:,j));
        if j==1, hold on; end
    end 
    hold off;
    shading interp;
    view(2);  
end


if shape==7
    r=hypot(X,Y);
    uex=(-6+7*r-r.^3)./(6*r);
    
    err=norm(uu(:)-uex(:),'inf')/norm(uex(:),'inf');
    disp(err);
    
    figure(2);
    for j=1:size(quads,1)
        surf(X(:,:,j), Y(:,:,j), uu(:,:,j)-uex(:,:,j));
        if j==1, hold on; end
    end 
    hold off;
    shading interp;
    view(2);  
end


if shape==8
    zet=X;
    rho=abs(Y);
    r1=hypot(rho,zet-z1);
    r2=hypot(rho,zet-z2);
    u1=m1*(R1./pi)^2*((r1<=R1).*(-1-R1*sinc(pi*r1/R1))+...
                      (r1> R1).*(-R1./(r1+(r1<=R1))));
    u2=m2*(R2./pi)^2*((r2<=R2).*(-1-R2*sinc(pi*r2/R2))+...
                      (r2> R2).*(-R2./(r2+(r2<=R2))));
    uex=u1+u2;
    err=norm(uu(:)-uex(:),'inf')/norm(uex(:),'inf');
    display(err);
    
    figure(2);
    for j=1:size(quads,1)
        surf(X(:,:,j), Y(:,:,j), uu(:,:,j)-uex(:,:,j));
        if j==1, hold on; end
    end 
    hold off;
    shading interp;
    view(2);
end

end

