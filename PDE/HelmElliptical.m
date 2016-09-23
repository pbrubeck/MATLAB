function [] = HelmElliptical(a, b, N, k)
f=sqrt(a^2-b^2);
[A1,B1,A2,B2,u,v]=chebLapEll(a,b,2*N-1,N);
A1=A1(2:end-1,2:end-1);
B1=B1(2:end-1,2:end-1);
C1=eye(2*N-3);
C2=-eye(N);

a=k^2;
lam=-2*a;
rmf=cos(pi/2*u(2:end-1)/u(end));
amf=cos(k*v(:));
[lam,a,X1,X2] = newton_mep(A1,B1,C1,A2,B2,C2,rmf,amf,lam,a);


target = [0 0];  		 % we are looking for the closest eigenvalue to the target
M1 = inv(full(A1-target(1)*B1-target(2)*C1)); % preconditioners
M2 = inv(full(A2-target(1)*B2-target(2)*C2)); 
solveM1 = @(x) M1*x;  % we can pass preconditioners as function handles
solveM2 = @(x) M2*x;

OPTS=[];							
OPTS.M1 = solveM1;    % preconditioners 
OPTS.M2 = solveM2;    
OPTS.minsize = 5;     
OPTS.maxsize = 10;
OPTS.maxsteps = 1000;   % maximum number of outer iterations
OPTS.extraction = 'mindist'; 
OPTS.reschange = 0; % 10^(-6); % we change to minimal residual when residual is less then epschange
OPTS.innersteps = 2;   % number of GMRES steps
OPTS.innertol = 1e-15;  % tolerance for the GMRES method
OPTS.target = target;
OPTS.delta = 5e-6;
OPTS.showinfo = 2; % set to 1 to see more information
OPTS.harmonic = 1;
OPTS.window = 0;
[lam,a,X1,X2] = twopareigs(A1,B1,C1,A2,B2,C2,k,OPTS);
q=-lam^2*f^2/4;

rmf=zeros(N,1);
rmf(2:end,:)=bsxfun(@times, conj(X1(N-1,:)), X1(1:N-1,:));
amf=bsxfun(@times, conj(X2(1,:)), X2);
rmf=rmf/max(rmf);
amf=amf/max(amf);
u=u(1:N);

display(a);
display(q);

figure(1); plot(u,rmf);
figure(2); plot(v,amf);

xx=f*cosh(u)*cos(v);
yy=f*sinh(u)*sin(v);
zz=rmf*amf.';

figure(3);
surfl(xx(:,[1:end 1]),yy(:,[1:end 1]),zz(:,[1:end 1]),'light'); 
shading interp;
axis off;
colormap(jet(256));
zrange=max(zz(:))-min(zz(:));
xrange=max(xx(:))-min(xx(:));
yrange=max(yy(:))-min(yy(:));
daspect([1 1 2*zrange/hypot(xrange,yrange)]);
view(0,90);
end