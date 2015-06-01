function [] = nbody( n )
% N-body problem simulation

% p is a 3-dimensional array whose page's columns are
% the position and velocity vectors for each particle.
p=zeros(3, 2, n);

a=1.0;
b=1.4;

% Initial state
for i=1:n
    th=2*pi*i/n;
    sine=sin(th);
    cosine=cos(th);
    p(:,1,i)=[a*cosine; b*sine; 0];
end
k=interact(p,n);
for i=1:n
    th=2*pi*i/n;
    sine=sin(th);
    cosine=cos(th);  
    
    g=k(:,2,i);
    A=[-a*cosine, -a*sine; -b*sine, b*cosine];
    w2=[1 0]*(A\g(1:2,:));
    p(:,2,i)=sqrt(w2)*[-a*sine; b*cosine; 0];
end

h=0.001;
frames=2000;
for i=1:frames
    % 4th order Classic Runge Kutta
    % p=RungeKutta(@(t,u) interact(u,n), p, 0, step, 1);
    k1=h*interact(p, n);
    k2=h*interact(p+k1/2, n);
    k3=h*interact(p+k2/2, n);
    k4=h*interact(p+k3, n);
    p=p+(k1+2*k2+2*k3+k4)/6;
    
    % Plot
    X=p(1,1,:);
    Y=p(2,1,:);
    Z=p(3,1,:);
    scatter3(X,Y,Z);
    drawnow;
end
end

% Gravitational interaction function
function k=interact(p, n)
k=zeros([3, 2, n]);
for i=1:n
    g=[0;0;0];
    for j=1:n
        if(i~=j)
            r=p(:,1,j)-p(:,1,i);
            r2=r'*r;
            g=g+r/sqrt(r2*r2*r2);
        end
    end
    k(:,1,i)=p(:,2,i);
    k(:,2,i)=g;
end
end