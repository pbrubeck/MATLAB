function [points,quads,curv] = bingrid(r1,r2,h)
a=1/6;
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
q0=[5,1,8,4; 6,2,5,1; 7,3,6,2; 8,4,7,3];
quad_inner=[q0+0;q0+4;q0+8;q0+12;  q0+20;q0+24;q0+28;q0+32];
% Middle
quad_mid=[38,37,17,18; 20,19,39,40; 18,37,19,40];
% Outer shells
quad_outer=[41,17,44,20; 42,38,41,17; 43,39,42,38; 44,20,43,39];
% Cores
quad_core=[1,2,4,3; 21,22,24,23];
% Gather all quads
quads=[quad_inner;quad_mid;quad_outer;quad_core];

% Radius of curvature
curv=zeros(size(quads));
curv(:,:)=inf;

% Inner right shells
curv(1 :4 ,1)=R_right(2);
curv(5 :8 ,2)=R_right(2);
curv(5 :8 ,1)=R_right(3);
curv(9 :12,2)=R_right(3);
curv(9 :12,1)=R_right(4);
curv(13:16,2)=R_right(4);
curv(15,1)=rc;
curv(13,1)=R1;

% Inner left shells
curv(17:20,1)=R_left(2);
curv(21:24,2)=R_left(2);
curv(21:24,1)=R_left(3);
curv(25:28,2)=R_left(3);
curv(25:28,1)=R_left(4);
curv(29:32,2)=R_left(4);
curv(29,1)=rc;
curv(31,1)=R1;

% Middle Patches
curv(33,1:2)=[R1, rc];
curv(34,1:2)=[R1, rc];
curv(35,3:4)=[-rc, rc];
curv(35,1:2)=[-rc, rc];

% Outer shell
curv(36:39,2)=R1;
curv(36:39,1)=R2;
end