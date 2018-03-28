function [points,quads,curv] = coregrid(r1)
R_right=[r1/4; r1/2; r1];
theta = 2*pi/4*(0.5:1:3.5)';

% Right shells
shell_right=exp(1i*theta)*R_right';
points=shell_right(:);

% Right and Left shells
q0=[5,1,8,4; 6,2,5,1; 7,3,6,2; 8,4,7,3];
quad_inner=[q0+0;q0+4];
% Core
quad_core=[1,2,4,3];
% Gather all quads
quads=[quad_inner;quad_core];

% Radius of curvature
curv=zeros(size(quads));
curv(:,:)=inf;

% Inner right shells
curv(1 :4 ,1)=R_right(2);
curv(5 :8 ,2)=R_right(2);
curv(5 :8 ,1)=R_right(3);
end

