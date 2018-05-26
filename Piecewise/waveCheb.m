function [] = waveCheb(N)
% Solves the wave equation using Chebyshev collocation
% spectral element method and Pade approximation of the matrix exponential.

% Subscript 0 : element / subdomain
% Subscript 1 : whole grid
% Subscript 2 : phase space (scalar field + spacetime momentum) [u,h,p]

% Problem parameters
vel=0;
mass=0;        % mass=0 yields wave equation
L=4*pi;        % Spatial domain [-L,L]
T=4*L;         % Time domain [0,T]
dt=0.1;
nsteps=ceil(T/dt);
dt=T/nsteps;

[D,x]=chebD(N);
x=L*x; x=x(end:-1:1);
D=D/L; D=D(end:-1:1,end:-1:1);
I=eye(N);

% Degrees of freedom
con=2:N-1; con=[];
rd=[con,N+1,3*N];     % Removed dofs
kd=setdiff(1:3*N,rd); % Kept dofs

% Constraint operator
B=zeros(length(rd),3*N);
B(1:end-2, 1:N)=D(con,:);
B(1:end-2, 2*N+1:3*N)=-I(con,:);
B(end-1,[N+1,2*N+1])=[1,1];
B(end,[2*N,3*N])=[1,-1];

% Schur complement basis, E1 satisfies homogenous BCs
E=eye(3*N,3*N);
E(rd,kd)=-B(:,kd);
E(rd,:)=B(:,rd)\E(rd,:);
E1=E(:,kd);
P=eye(3*N,3*N);
P1=P(:,kd);

% Time propagator
maskM1=[0,-1,0; mass,0,0; 0,0,0];
maskK1=[0,0,0; 0,0,-1; 0,-1,0];
A=P1'*(kron(maskM1,I)+kron(maskK1,D))*E1;
U=expm(A*dt);

figure(4);
imagesc(log(abs(A))); colormap(gray(256)); colorbar();

[Lambda]=eig(U,'vector');
figure(5);
plot(Lambda,'.');

% Initial condition
u=-1/2*sin((-vel*x-abs(x))/(1-vel^2));
h= 1/2*cos((-vel*x-abs(x))/(1-vel^2)).*(1+vel*sign(x))/(1-vel^2);
p= 1/2*cos((-vel*x-abs(x))/(1-vel^2)).*(vel+1*sign(x))/(1-vel^2);
q=[u; h; p];

% Plot initial condition
set(0,'defaultfigureposition',[0,0,1080,720]);
figure(1); h1=plot(x,u,'b','Linewidth',2); xlim([-L,L]); ylim([-0.5,0.5]); title('Scalar wave');
figure(2); h2=plot(x,h,'b','Linewidth',2); xlim([-L,L]); ylim([-0.5,0.5]); title('Enthalpy');
figure(3); h3=plot(x,p,'b','Linewidth',2); xlim([-L,L]); ylim([-0.5,0.5]); title('Momentum');

% Jump Force
xi=0;
Jumps=zeros(N+1,3);
Ju=zeros(N,1); Ju(2:4:end)=1;  Ju(4:4:end)=-1;
Jh=zeros(N,1); Jh(2:4:end)=1i; Jh(4:4:end)=-1i;
Jp=zeros(N,1); Jp(1:4:end)=1;  Jp(3:4:end)=-1;
Jumps(1:N,:)=[Ju,Jh,Jp];

JK=zeros(N,3);
JK(:,2)=Jumps(1:N,3);
JK(:,3)=Jumps(1:N,2);

id1=1:N/2;
id2=N/2+1:N;
[s1,s2]=piecewiseLagrange(x/L,xi/L,JK);
f=zeros(N,3);
f(id1,:)=-D(id1,:)*s1;
f(id2,:)=-D(id2,:)*s2;
f=f(:);

%Steady state
b=-(A+1i*eye(size(A,1)))\(P1'*f);

% Exact solution
uex=@(t) -1/2*sin((t-vel*x-abs(x-vel*t))/(1-vel^2));

% Create .gif file
filename='waveCheb.gif';
im=frame2im(getframe(1));
[imind,cm]=rgb2ind(im,256);
imwrite(imind,cm,filename,'gif','DelayTime',0,'Loopcount',inf);

% Time-stepping
t=0;
for i=1:nsteps
    t=t+dt;
    q=E1*(U*real(q(kd)-b)+real(b*exp(-1i*dt)));
    b=b*exp(-1i*dt);
    
    set(h1,{'YData'},{q(1:N)});
    set(h2,{'YData'},{q(N+1:2*N)});
    set(h3,{'YData'},{q(2*N+1:3*N)});
    drawnow;
    
    im=frame2im(getframe(1));
    [imind,cm]=rgb2ind(im,256);
    imwrite(imind,cm,filename,'gif','DelayTime',0.001,'WriteMode','append');
end

end