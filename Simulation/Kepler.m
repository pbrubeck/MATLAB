function [t]=Kepler(t)
% Orbital parameters for the Earth and Mars

% Semi-latus rectum (AU)
L1=0.999734818;
L2=1.51025679;

% Eccentricity
ec1=0.016710220;
ec2=0.093412330;

% Angular velocity (rad/year)
w1=2*pi;
w2=2*pi/1.8808476;

% True anomaly 01/03/15 (rad)
phi1=6.265;
phi2=0.2421;

% Orthonormal Basis
% [Node, Inclination, Longitude]
Rot1=EulerRotation(-0.1965, 9E-07, 1.7967674);
Rot2=EulerRotation(0.8653, 0.0322992, 5.8650191);


k=6.2832;
f1=@(y,t) k*sqrt(((1-ec1*cos(y))/L1)^3);
f2=@(y,t) k*sqrt(((1-ec2*cos(y))/L2)^3);

n=200;
h=1.0000000/n;
th1=RungeKuta(f1, phi1, 0, h, n);
h=1.8808476/n;
th2=RungeKuta(f2, phi2, 0, h, n);

hold on
plot(linspace(0,1,n), th1);
plot(linspace(0,1.8808476,n), th2);
hold off

ii=0;
err=0;
while(ii<200 && (ii==0 || abs(err)>1E-15))
    [r1, v1, a1]=orbit(L1, ec1, phi1, w1, t, Rot1);
    [r2, v2, a2]=orbit(L2, ec2, phi2, w2, t, Rot2);
    s0=r1-r2;
    s1=v1-v2;
    s2=a1-a2;
    
    f1=2*s0'*s1;
    f2=2*(s1'*s1+s0'*s2);
    fprintf('i=%d \t t=%f \t |s|=%f \t f''(t)=%e \n', ii, t, norm(s0), f1);
    err=f1/f2;
    t=t-err;
    ii=ii+1;
end
end

function [pos, vel, acc]=orbit(L, ec, phi, w, t, Rot)
angle=w*t+phi;
sine=sin(angle);
cosine=cos(angle);

r0=L/(1+ec*cosine);
r1=w*ec*sine*r0^2/L;
r2=w*ec*r0*(w*r0*cosine+2*r1*sine)/L;

er=Rot*[ cosine; sine; 0];
et=Rot*[-sine; cosine; 0];

pos=r0*er;
vel=r1*er+w*r0*et;
acc=(r2-w*w*r0)*er+2*w*r1*et;
end