function [] = schrodExpTiles( m )
% Time-dependent Schrodinger equation with domain decomposition method 
% using Schur complement
n=m;

% Timestep
dt=0.01;

% L-shaped membrane
adjx=[2 1];
adjy=[3 1];

% UIUC block-I
adjx=[3 4; 4 5; 6 7; 7 8];
adjy=[4 2; 2 1; 1 7];

% Topology
[topo,net,RL,TB]=ddtopo(adjx,adjy);
pos=ddpatches(topo);
x0=real(pos);
y0=imag(pos);

% Degrees of freedom
rd1=[1,m];
rd2=[1,n];
kd1=2:m-1;
kd2=2:n-1;

E1=eye(m);
E2=eye(n);

% Differential operators
[Dx,x]=chebD(m);
[Dy,y]=chebD(n);
H1=-1/2*Dx*Dx;
H2=-1/2*Dy*Dy;

A1=myexp(-dt/2i*H1);
A2=myexp(-dt/2i*H2);

% Constraint operators
C1=E1(rd1,:);
C2=E2(rd2,:);

% Left Propagator Schur complement
[S,~,V1,V2,Q]=ddschurprod(adjx,adjy,A1,A2,Dx,Dy,C1,C2);
[Lschur, Uschur, bschur]=lu(S,'vector');

figure(2);
imagesc(log(abs(S)));
colormap(gray(256)); colorbar;
drawnow;

% Propagation in a square [-1,1]^2
W1=inv(V1(kd1,kd1));
W2=inv(V2(kd2,kd2));
lhsprop=@(uu) V1(kd1,kd1)*(Q(kd1,kd2).*(W1*uu*W2.'))*V2(kd2,kd2).';

% Right propagator
rhspropfull=@(uu) conj(A1)*uu*conj(A2).';

% Poisson solver
function [uuu]=propagate(uuu)
    F=zeros(m-2,n-2,size(uuu,3));
    for j=1:size(uuu,3)
        pp=rhspropfull(uuu(:,:,j));
        F(:,:,j)=pp(kd1,kd2);
    end
    
    v=cell([size(net,1),1]);
    for j=1:size(net,1)
        v{j}=lhsprop(F(:,:,j));
    end
    
    rhs=zeros(m-2, size(RL,1)+size(TB,1));
    for j=1:size(RL,1)
        rhs(:,RL(j,1))=-(Dx(rd1(2),kd1)*v{adjx(j,1)}-Dx(rd1(1),kd1)*v{adjx(j,2)});
    end
    for j=1:size(TB,1)
        rhs(:,TB(j,1))=-(v{adjy(j,1)}*Dy(rd2(2),kd2)'-v{adjy(j,2)}*Dy(rd2(1),kd2)');
    end
    rhs=rhs(:);

    % Solve for boundary nodes
    b=Uschur\(Lschur\rhs(bschur));
    b=reshape(b, m-2, []);
    b=[b, zeros(m-2,1)];
    
    % Solve for interior nodes with the given BCs
    for j=1:size(uuu,3)
        uuu(rd1,kd2,j)=b(:,net(j,1:2)).';
        uuu(kd1,rd2,j)=b(:,net(j,3:4));
        uuu(kd1,kd2,j)=lhsprop(F(:,:,j)-A1(kd1,rd1)*b(:,net(j,1:2)).'-b(:,net(j,3:4))*A2(kd2,rd2).');        
    end
end

[xx,yy]=ndgrid(x,y);
uuu=zeros(m,n,size(net,1));
sig=0.2;
zg=[-1i; 3i];

for i=1:size(uuu,3)
    for k=1:size(zg,1)
        uuu(:,:,i)=uuu(:,:,i)+exp(-((xx+x0(i)-real(zg(k))).^2+(yy+y0(i)-imag(zg(k))).^2)/(2*sig^2));
    end
end

figure(1);
h=cell(size(uuu,3),1);
for i=1:size(uuu,3)
    h{i}=surf(xx+x0(i), yy+y0(i), real(uuu(:,:,i)));
    if i==1, hold on; end;
end
hold off;
colormap(jet(256));  view(2); 
shading interp;
zlim([-1,1]);
axis square manual;
camlight;

function timeEvol(x,y)
    uuu=propagate(uuu);
    for p=1:size(uuu,3)
        set(h{p}, 'Zdata', real(uuu(:,:,p)));
    end
    drawnow;
end
fps=20;
tf=10;
nplots=ceil(tf/dt);
t=timer('StartDelay', 0, 'Period', round(1/fps,3), 'TasksToExecute', nplots, 'ExecutionMode','fixedRate');
t.TimerFcn = @timeEvol;
start(t);
end