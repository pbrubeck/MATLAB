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

h=0.001 ;
frames=1000;
for i=1:frames
    % Runge Kutta Felhberg 45
    
    k1=h*interact(p, n);
    k2=h*interact(p+1/4*k1, n);
    k3=h*interact(p+3/32*k1+9/32*k2, n);
    k4=h*interact(p+1932/2197*k1-7200/2197*k2+7296/2197*k3, n);
    k5=h*interact(p+439/216*k1-8*k2+3680/513*k3-845/4104*k4, n);
    k6=h*interact(p-8/27*k1+2*k2-3544/2565*k3+1859/4104*k4-11/40*k5, n);
    
    R=errNorm(1/360*k1-128/4275*k3-2197/75240*k4+1/50*k5+2/55*k6, n)/h;
    q=(1E-5/(2*R))^(1/4);
    if R<1E-5
       
    end
    p=p+25/216*k1+1408/2565*k3+2197/4104*k4-1/5*k5;

    % Plot
    X=p(1,1,:);
    Y=p(2,1,:);
    Z=p(3,1,:);
    figure(1);
    scatter3(X,Y,Z, 'filled');
    axis([-1 1 -1 1 -1 1]);

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

% Error norm
function err=errNorm(k, n)
err=0;
for i=1:n
    x=k(:,1,i);
    v=k(:,2,i);
    err=err+x'*x+v'*v;
end
err=sqrt(err);
end
