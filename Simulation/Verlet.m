function [] = Verlet( n )
% N-body problem simulation

% p is a 3-dimensional array whose pages are the position 
% and velocity vectors for each particle.
p=zeros(3, n, 2);

a=1.0;
b=1.0;
h=0.001; h2=h*h;

% Initial state
for i=1:n
    th=2*pi*i/n;
    sine=sin(th);
    cosine=cos(th);
    p(:,i,2)=[a*cosine; b*sine; 0];
end
k=interact(p(:,:,2),n);
for i=1:n
    th=2*pi*i/n;
    sine=sin(th);
    cosine=cos(th);  
    A=[-a*cosine, -a*sine; -b*sine, b*cosine];
    w2=[1 0]*(A\k(1:2,i));
    p(:,i,1)=p(:,i,2)+h*sqrt(w2)*[-a*sine; b*cosine; 0];
end

frames=2000;
plot=scatter3(p(1,:,1), p(2,:,1), p(3,:,1));
for i=1:frames
    temp=p(:,:,1);
    delta=p(:,:,1)-p(:,:,2);
    k1=delta+h2*interact(temp, n);
    k2=delta+h2*interact(temp+k1/2, n);
    k3=delta+h2*interact(temp+k2/2, n);
    k4=delta+h2*interact(temp+k3, n);
    p(:,:,1)=temp+(k1+2*k2+2*k3+k4)/6;
    p(:,:,2)=temp;
    
    % Plot
    plot.XData=p(1,:,1);
    plot.YData=p(2,:,1);
    plot.ZData=p(3,:,1);
    drawnow;
end
end

% Gravitational interaction function
function k=interact(p, n)
k=zeros(3, n);
for i=1:n
    g=[0;0;0];
    for j=1:n
        if(i~=j)
            r=p(:,j)-p(:,i);
            r2=r'*r;
            g=g+r/sqrt(r2*r2*r2);
        end
    end
    k(:,i)=g;
end
end