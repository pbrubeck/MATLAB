r1 = 1;
r2 = 2.5;
h = 5;
k = 1.5;

y0 = 2*max(r1,r2);
R0 = hypot(y0,h);
b = k*hypot(y0,h);
a = k*(h+hypot(y0,h));
theta = zeros(8,1);
theta(1:2:8) = (0:3)*pi/2;
theta(2) = atan2(y0,h);
theta(4) = atan2(y0,-h);
theta(6) = atan2(-y0,-h);
theta(8) = atan2(-y0,h);

phi = pi/6*(0:11)';

rightC = h+r1*exp(1i*theta); % inner circle right side
leftC = -h+r2*exp(1i*theta); % inner circle left side
R =  h+R0*exp(1i*theta); % mid circle right
S = -h+R0*exp(1i*theta); % mid circle left
outerC = [R(1:4); S(3:7); R(6:8)]; % outermost circle
outerE = a.*exp(1i*phi);

points = [rightC;leftC;outerC;outerE;0];
figure(1); clf;
% plot(points,'*');
% axis equal;

nquad = 28;
quad = zeros(nquad,4);
quad(:,4) = 1:28;
quad(:,3) = quad([2:8,1,10:16,9,18:28,17],4);
quad(:,1) = [18:20,41,26:28,17,20:26,41,30:40,29];
quad(:,2) = quad([8,1:7,16,9:15,28,17:27],1);

% fence(i,:) refers gives information on interface i.
% key: [East, West, North, South] = [1, 2, 3, 4].
% i.e. fence(i,:)= [3,8,4,1].
% interface i is on North (3) of Quad 8, and South (4) of Quad 1.
% fence(i,1) is the side of quad 1 found on interface i.
% fence(i,2) is quad 1.
% fence(i,3) is the side of quad 2 found on interface i.
% fence(i,4) is quad 2.
nfence = 42;
fence = zeros(nfence,4);
fence(1:28,1) = 4;
fence(1:28,2) = 1:28;
fence(1:28,3) = 3;
fence(1:28,4)= [8,1:7,16,9:15,28,17:27]; %north on radial boundries
fence(29:40,1) = 2;
fence(29:40,2) = 17:28;
fence(29:40,3) = 1;
fence(29:40,4) = [1:3,10:15,6:8];
fence(41,:) = [1,4,1,9];
fence(42,:) = [1,5,1,16];


X=real(points(quad(:,[2,1,3,4])'));
Y=imag(points(quad(:,[2,1,3,4])'));
colormap(jet(size(quad,1)));
C = 1:size(quad,1);
patch(X,Y,C)
axis equal;
grid on;