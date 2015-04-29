function [] = nbody( n )
% N-body problem simulation

% p is a 3-dimensional array whose page's columns are
% the position and velocity vectors for each particle.
p=zeros([3, 2, n]);

% Initial state
for i=1:n
    th=2*pi*i/n;
    p(:,1,i)=[cos(th);sin(th);0];
end
k=interact(p,n);
for i=1:n
    r=p(:,1,i);
    g=k(:,2,i);
    p(:,2,i)=sqrt(-g'*r)*[0 -1 0;1 0 0;0 0 0]*r;
end

h=0.001;
frames=1000;
for i=1:frames
    % 4th order Classic Runge Kutta
    % p=RungeKutta(@(t,u) interact(u,n), p, 0, step, 1);
    k1=interact(p, n);
    k2=interact(p+h/2*k1, n);
    k3=interact(p+h/2*k2, n);
    k4=interact(p+h*k3, n);
    p=p+h*(k1+2*k2+2*k3+k4)/6;
    
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
            g=g+r/sqrt((r'*r)^3);
        end
    end
    k(:,1,i)=p(:,2,i);
    k(:,2,i)=g;
end
end