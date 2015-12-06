function [] = gustavus(x0,y0,vx0,vy0,b,h,n)
g=[0; -9.81];
pos=zeros(2,n);
vel=zeros(2,n);
pos(:,1)=[x0; y0];
vel(:,1)=[vx0; vy0];
for i=1:n
    acc=b*vel(:,i)+g;
    vel(:,i+1) = vel(:,i)+h*acc;
    pos(:,i+1) = pos(:,i)+h*(vel(:,i)+0.5*h*acc);
end
plot(pos(1,:), pos(2,:));
disp(pos(:,end));
end