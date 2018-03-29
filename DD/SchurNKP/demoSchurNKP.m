function [lam] = demoSchurNKP(shape, ref, m, k)
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
        case 3
            [z0,quads,curv]=bingrid(1,1,3);
    end
end

% Mesh refinement
for j=1:ref
    [~, adj, ~, ~, bnd] = meshtopo(quads);
    [z0,quads,curv]=quadmeshrefine(z0,quads,curv,adj,bnd);
end

if nargin>3
    lam = meshSchurNKP(z0, quads, curv, m, k);
else
    lam = meshSchurNKP(z0, quads, curv, m);
end

% delete(findall(gcf,'Type','light'));
% shading faceted;
end

