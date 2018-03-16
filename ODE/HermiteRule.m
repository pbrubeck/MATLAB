function [] = HermiteRule(N)
% Time integration with generalized Hermite rule.

% du/dt = f(t)
f=@(t) sin(t);
% Initial value
u0=0;
% Exact solution
uex=@(t) 1-cos(t);

% Time interval
t0=0;
tf=30*pi;

maxorder = 6;
err = zeros(length(N),maxorder);
tms = zeros(length(N),maxorder);

% Timestep
for i=1:length(N)
    h=(tf-t0)/(N(i)-1);
    t=(t0:h:tf)';
    
    tic;
    u2 = herm2(f,t,u0);
    tms(i,1)=toc();
    err(i,1)=norm(u2-uex(t))/numel(t);
    
    tic;
    u4 = herm4(f,t,u0);
    tms(i,2)=toc();
    err(i,2)=norm(u4-uex(t))/numel(t);

    tic;
    u6 = herm6(f,t,u0);
    tms(i,3)=toc();
    err(i,3)=norm(u6-uex(t))/numel(t);

    tic;
    u8 = herm8(f,t,u0);
    tms(i,4)=toc();
    err(i,4)=norm(u8-uex(t))/numel(t);
    
    tic;
    u10 = herm10(f,t,u0);
    tms(i,5)=toc();
    err(i,5)=norm(u10-uex(t))/numel(t);
    
    tic;
    u12 = herm12(f,t,u0);
    tms(i,6)=toc();
    err(i,6)=norm(u12-uex(t))/numel(t);
    
    tic;
    u14 = herm14(f,t,u0);
    tms(i,7)=toc();
    err(i,7)=norm(u14-uex(t))/numel(t);
end

figure(1);
loglog(N,err);

figure(2);
loglog(err,tms,'-*');
end


function u=herm2(f,t,u0)
% 2nd order Hermite rule
difforder = 1;
t1 = ainit(t(1:end-1), difforder);
t2 = ainit(t(2:end),   difforder);
f1 = f(t1);
f2 = f(t2);
du=zeros(size(t));
h = diff(t);
du(2:end) = (h.^1/2).*(f1{0}+f2{0});
u=u0+cumsum(du);
end

function u=herm4(f,t,u0)
% 4th order Hermite rule
difforder = 1;
t1 = ainit(t(1:end-1), difforder);
t2 = ainit(t(2:end),   difforder);
f1 = f(t1);
f2 = f(t2);
du=zeros(size(t));
h = diff(t);
du(2:end) = (h.^1/2 ).*(f1{0}+f2{0})+...
            (h.^2/12).*(f1{1}-f2{1});
u=u0+cumsum(du);
end

function u=herm6(f,t,u0)
% 6th order Hermite rule
difforder = 2;
t1 = ainit(t(1:end-1), difforder);
t2 = ainit(t(2:end),   difforder);
f1 = f(t1);
f2 = f(t2);
du=zeros(size(t));
h = diff(t);
du(2:end) = (h.^1/2  ).*(f1{0}+f2{0})+...
            (h.^2/10 ).*(f1{1}-f2{1})+...
            (h.^3/120).*(f1{2}+f2{2});
u=u0+cumsum(du);
end

function u=herm8(f,t,u0)
% 8th order Hermite rule
difforder = 3;
t1 = ainit(t(1:end-1), difforder);
t2 = ainit(t(2:end),   difforder);
f1 = f(t1);
f2 = f(t2);
du=zeros(size(t));
h = diff(t);
du(2:end) = (h.^1/2   ).*(f1{0}+f2{0})+...
            (h.^2*3/28).*(f1{1}-f2{1})+...
            (h.^3/84  ).*(f1{2}+f2{2})+...
            (h.^4/1680).*(f1{3}-f2{3});
u=u0+cumsum(du);
end

function u=herm10(f,t,u0)
% 10th order Hermite rule
difforder = 4;
t1 = ainit(t(1:end-1), difforder);
t2 = ainit(t(2:end),   difforder);
f1 = f(t1);
f2 = f(t2);
du=zeros(size(t));
h = diff(t);
du(2:end) = (h.^1/2    ).*(f1{0}+f2{0})+...
            (h.^2/9    ).*(f1{1}-f2{1})+...
            (h.^3/72   ).*(f1{2}+f2{2})+...
            (h.^4/1008 ).*(f1{3}-f2{3})+...
            (h.^5/30240).*(f1{4}+f2{4});
u=u0+cumsum(du);
end

function u=herm12(f,t,u0)
% 12th order Hermite rule
difforder = 5;
t1 = ainit(t(1:end-1), difforder);
t2 = ainit(t(2:end),   difforder);
f1 = f(t1);
f2 = f(t2);
du=zeros(size(t));
h = diff(t);
du(2:end) = (h.^1*1/2     ).*(f1{0}+f2{0})+...
            (h.^2*5/44    ).*(f1{1}-f2{1})+...
            (h.^3*1/66    ).*(f1{2}+f2{2})+...
            (h.^4*1/792   ).*(f1{3}-f2{3})+...
            (h.^5*1/15840 ).*(f1{4}+f2{4})+...
            (h.^6*1/665280).*(f1{5}-f2{5});
u=u0+cumsum(du);
end

function u=herm14(f,t,u0)
% 12th order Hermite rule
difforder = 6;
t1 = ainit(t(1:end-1), difforder);
t2 = ainit(t(2:end),   difforder);
f1 = f(t1);
f2 = f(t2);
du=zeros(size(t));
h = diff(t);
du(2:end) = (h.^1*1/2       ).*(f1{0}+f2{0})+...
            (h.^2*3/26      ).*(f1{1}-f2{1})+...
            (h.^3*5/312     ).*(f1{2}+f2{2})+...
            (h.^4*5/3432    ).*(f1{3}-f2{3})+...
            (h.^5*1/11440   ).*(f1{4}+f2{4})+...
            (h.^6*1/308880  ).*(f1{5}-f2{5})+...
            (h.^7*1/17297280).*(f1{6}+f2{6});
u=u0+cumsum(du);
end