function []=collision()
L=1;

a= [1; 3; 0];
da=[1; 6; 1];

b= [1; 1; 1];
db=[1; 0; -1];

s=a-b;
u=da-db;
uu=u'*u;
su=s'*u;
ss=s'*s;
t=-(su+sqrt(su^2-uu*(ss-L^2)))/uu;
r=s+u*t;
dv=r*(u'*r)/(r'*r);
disp(r);

a=a+dv*t;
da=da+dv;
b=b-dv*t;
db=db-dv;

disp([a b]);
end