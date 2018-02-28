r1 = 2.5;
r2 = 2.5;
h = 5;
k = 1.5;

y0 = 2*max(r1,r2);
R0 = hypot(y0,h);
R1 = 2*h+hypot(y0,h);
theta = zeros(8,1);
theta(1:2:8) = (0:3)*pi/2;
theta(2) = atan2( y0, h);
theta(4) = atan2( y0,-h);
theta(6) = atan2(-y0,-h);
theta(8) = atan2(-y0, h);
phi = pi/6*(0:11)';

rightC = h+r1*exp(1i*theta);     % right inner circle
leftC = -h+r2*exp(1i*theta);     % left inner circle
R =  h+R0*exp(1i*theta);         % mid circle right
S = -h+R0*exp(1i*theta);         % mid circle left
midC = [R(1:4); S(3:7); R(6:8)]; % mid circles
outC = R1.*exp(1i*phi);          % outer circle
points = [rightC;leftC;midC;outC;0;h;-h];

nquad = 36;
quad = zeros(nquad,4);
quad(:,4) = [1:28, 42,42,42,42, 43,43,43,43];                 % 11 SW
quad(:,1) = [18:20,41, 26:28,17, 20:26,41, 30:40,29, 2:2:16]; % 00 NE
quad(:,2) = [quad([2:8,1, 10:16,9, 18:28,17],4)',    1:2:16]; % 01 NW
quad(:,3) = [quad([8,1:7, 16,9:15, 28,17:27],1)',    3:2:8,1, 11:2:16,9]; % 10 SE

% edge(i,:) refers gives information on interface i.
% key: [East, West, North, South] = [1, 2, 3, 4].
% i.e. fence(i,:)= [3,8,4,1].
% interface i is on North (3) of Quad 8, and South (4) of Quad 1.
% edge(i,1) is the side of quad 1 found on interface i.
% edge(i,2) is quad 1.
% edge(i,3) is the side of quad 2 found on interface i.
% edge(i,4) is quad 2.

nedge = 42;
edge = zeros(nedge,4);
edge(1:28, 1) = 4;
edge(1:28, 2) = 1:28;
edge(1:28, 3) = 3;
edge(1:28, 4) = [8,1:7, 16,9:15, 28,17:27];
edge(29:40,1) = 2;
edge(29:40,2) = 17:28;
edge(29:40,3) = 1;
edge(29:40,4) = [1:3, 10:15, 6:8];
edge(41,:) = [1,4,1,9];
edge(42,:) = [1,5,1,16];

% Radius of curvature
curv=zeros(nquad,4);
curv(:,:)=inf;
curv([1:3,6:8],1)=R0;
curv(1:8,  2)=r1;
curv(10:15,1)=R0;
curv(9:16, 2)=r2;
curv(17:28,1)=R1;
curv(17:28,2)=R0;
curv(29:32,1)=-r1;
curv(29:32,3)= r1;
curv(33:36,1)=-r2;
curv(33:36,3)= r2;

% Sample gird
N=16;
x=linspace(-1,1,N)';
y=linspace(-1,1,N)';
[xx,yy]=ndgrid(x,y);

f0=@(r) exp(-r.^2/2).*(1-r.^2);
f=@(z) f0(abs(z-h)/r1)+f0(abs(z+h)/r2);

figure(1);  
for k=1:nquad
    F = curvedquad(points(quad(k,:)), curv(k,:));
    [jac,G11,G12,G22] = diffgeom(F,x,y); 
    zz = F(xx,yy);
    surf(real(zz), imag(zz), jac*0+k); hold on;
end
colormap(prism(nquad));
hold off;
axis square;
% shading interp
% camlight;
view(2);