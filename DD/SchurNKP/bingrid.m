r1 = 1;
r2 = 1;
h = 1.8;

a=1/4;
r0=a*h+(1-a)*max(r1,r2);
rc=h/sqrt(2);
R1=sqrt(2)*h+r0;
R2=2*R1;

R_right=[r1/4; r1/2; r1; r0; rc];
R_left =[r2/4; r2/2; r2; r0; rc];
theta = 2*pi/4*(0.5:1:3.5)';

% Right shells
shell_right=h+exp(1i*theta)*R_right';
shell_right([17,20])=R1*exp(1i*theta([1,4]));
% Left shells
shell_left=-h+exp(1i*theta)*R_left';
shell_left([18,19])=R1*exp(1i*theta([2,3]));
% Outer shell
shell_outer=R2/R1*[shell_right(end-3);shell_left([end-2;end-1]);shell_right(end)];
% Gather all points
points=[shell_right(:); shell_left(:); shell_outer(:)];

% Right and Left shells
q0=[5,8,1,4; 6,5,2,1; 7,6,3,2; 8,7,4,3];
quad_inner=[q0+0;q0+4;q0+8;q0+12;  q0+20;q0+24;q0+28;q0+32];
% Middle
quad_mid=[38,17,37,18; 20,39,19,40; 18,19,37,40];
% Outer shells
quad_outer=[41,44,17,20; 42,41,38,17; 43,42,39,38; 44,43,20,39];
% Cores
quad_core=[1,4,2,3; 21,24,22,23];
% Gather all quads
quad=[quad_inner;quad_mid;quad_outer;quad_core];
nquad=size(quad,1);

% Radius of curvature
curv=zeros(nquad,4);
curv(:,:)=inf;
% Inner right shells
curv(1 :4 ,3)=R_right(2);
curv(5 :8 ,4)=R_right(2);
curv(5 :8 ,3)=R_right(3);
curv(9 :12,4)=R_right(3);
curv(9 :12,3)=R_right(4);
curv(13:16,4)=R_right(4);
curv(13,3)=R1;
curv(15,3)=rc;

% Inner left shells
curv(17:20,3)=R_left(2);
curv(21:24,4)=R_left(2);
curv(21:24,3)=R_left(3);
curv(25:28,4)=R_left(3);
curv(25:28,3)=R_left(4);
curv(29:32,4)=R_left(4);
curv(31,3)=R1;
curv(29,3)=rc;

% Middle Patches
curv(33,3:4)=[R1, rc];
curv(34,3:4)=[R1, rc];
curv(35,1:2)=[-rc, rc];
curv(35,3:4)=[-rc, rc];

% Outer shell
curv(36:39,4)=R1;
curv(36:39,3)=R2;

% Sample gird
N=12;
x=linspace(-1,1,N)';
y=linspace(-1,1,N)';
[xx,yy]=ndgrid(x,y);

f0=@(r) exp(-r.^2/2).*(1-r.^2);
f=@(z) f0(abs(z-h)/r1)+f0(abs(z+h)/r2);

figure(1); clf;
for k=1:nquad
    F = curvedquad(points(quad(k,:)), curv(k,:));
    [jac,G11,G12,G22] = diffgeom(F,x,y); 
    zz = F(xx,yy);
    dd = abs(zz-h)-abs(zz+h);
    surf(real(zz), imag(zz), dd); hold on;
end
colormap(jet(256));
hold off;
axis square;
%shading interp
%camlight;
view(2);
