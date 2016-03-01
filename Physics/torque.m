function [] = torque(  )
origin=zeros(3,4);
r=eye(3);

figure(1); 
h=quiver3(origin(1,:),origin(2,:),origin(3,:),origin(1,:),origin(2,:),origin(3,:));
set(h,'Color', [1 0 0]);
axis equal; xlim([-1,1]); ylim([-1,1]); zlim([-1,1]); 


n=1024;
dt=2*pi/n;
for i=0:n-1
    t=i/(n-1);
    omega=-1;
    dtheta=omega*dt;
    
    phi=4*pi*t;
    theta=pi/3;
    rot=[cos(phi)*sin(theta); sin(phi)*sin(theta); cos(theta)];
    
    dq=[cos(dtheta); sin(dtheta)*rot/norm(rot)];
    r=quatrotate(dq', r);
    
    v=[r, rot];
    set(h,'UData',v(1,:),'VData',v(2,:),'WData',v(3,:));
    drawnow();
end
end