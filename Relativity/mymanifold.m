function [G11,G21,G12,G22,K11,K21,K12,K22,vol] = mymanifold(u,v,iflag)
m=length(u);
n=length(v);

% initialize audi grid
[u,v] = ndgrid(u,v);
[u,v] = ainit(u,v,2);

% define surface parametrization
x = u;
y = v;
z = exp(-10*(u.^2+v.^2));
S = [x;y;z];

% Differential geometry
J = ajac(S);                              % Jacobian
N = cross(J(:,1),J(:,2));                 % normal vector
N = N/norm(N);                            % normalize
G = J'*J;                                 % first fundamental form
K = [N'*adiff(S,2,0) N'*adiff(S,1,1);     % second fundamental form
     N'*adiff(S,1,1) N'*adiff(S,0,2)];

detG=det(G);

% Evaluation
if strcmp(iflag, 'inv')
    invG=inv(G);
    gg=invG{0};
else
    gg=G{0};
end
kk=K{0};

G11=reshape(gg(1,1,:), m,n);
G21=reshape(gg(2,1,:), m,n);
G12=reshape(gg(1,2,:), m,n);
G22=reshape(gg(2,2,:), m,n);
K11=reshape(kk(1,1,:), m,n);
K21=reshape(kk(2,1,:), m,n);
K12=reshape(kk(1,2,:), m,n);
K22=reshape(kk(2,2,:), m,n);
vol=reshape(sqrt(detG{0}), m,n);
end
