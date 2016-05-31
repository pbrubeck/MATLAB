function [] = EWE(N,m)
% Ellipsoidal Wave (Helmholtz) Equation
a=2;
b=sqrt(2);
r=4; % x3 boundary

[D,x]=chebD(N);

x1=b*x; 
D1=1/b*D;
L1=diag((x1.^2-a^2).*(x1.^2-b^2))*D1^2+diag(x1.*(2*x1.^2-(a^2+b^2)))*D1;

x2=b+(a-b)/2*(x+1);
D2=2/(a-b)*D;
L2=diag((x2.^2-a^2).*(x2.^2-b^2))*D2^2+diag(x2.*(2*x2.^2-(a^2+b^2)))*D2;

x3=a+(r-a)/2*(x+1);
D3=2/(r-a)*D;
L3=diag((x3.^2-a^2).*(x3.^2-b^2))*D3^2+diag(x3.*(2*x3.^2-(a^2+b^2)))*D3;

A={L1(2:end-1,2:end-1), L2(2:end-1,2:end-1), L3(2:end-1,2:end-1)};
B={diag(x1(2:end-1).^4), diag(x2(2:end-1).^4), diag(x3(2:end-1).^4)};
C={diag(x1(2:end-1).^2), diag(x2(2:end-1).^2), diag(x3(2:end-1).^2)};
D={eye(N-2), eye(N-2), eye(N-2)};

u1=cos(m*pi/b*x1);
u2=cos(m*pi/(a-b)*(x2-b));
u3=cos(m*pi/(r-a)*(x3-a));
u={u1(2:end-1), u2(2:end-1), u3(2:end-1)};
lam=[-m^2;-m^2;-m^2];

[u, lam]=mep(A,B,C,D,u,lam);
u1=[0; u{1}; 0];
u2=[0; u{2}; 0];
u3=[0; u{3}; 0];

figure(1); plot(x1, u1);
figure(2); plot(x2, u2);
figure(3); plot(x3, u3);
disp(lam);
end


function [u,lam]=mep(A,B,C,D,u,lam)
% A, B, C are arrays of matrices
n1=size(A{1},1);
n2=size(A{2},1);
n3=size(A{3},1);
n=n1+n2+n3;
id1=1:n1;
id2=1+n1:n1+n2;
id3=1+n1+n2:n;

x=[normc(u{1}); normc(u{2}); normc(u{3}); lam(:)];

P=A;
dx=ones(size(x));
J=zeros(size(x,1));

while(max(abs(dx))>1E-10)
    u1=x(id1);
    u2=x(id2);
    u3=x(id3);
    lam=x(end-2:end);
    for i=1:3
        P{i}=A{i}-lam(1)*B{i}-lam(2)*C{i}-lam(3)*D{i};
    end
    J(id1, id1)=P{1};
    J(id2, id2)=P{2};
    J(id3, id3)=P{3};
    J(end-2, id1)=2*u1';
    J(end-1, id2)=2*u2';
    J(end-0, id3)=2*u3';
    J(1:n, end-2:end)=[
        -B{1}*u1, -C{1}*u1, -D{1}*u1;
        -B{2}*u2, -C{2}*u2, -D{2}*u2;
        -B{3}*u3, -C{3}*u3, -D{3}*u3];
    f=[P{1}*u1; P{2}*u2; P{3}*u3; u1'*u1-1; u2'*u2-1; u3'*u3-1];
    dx=J\f;
    x=x-dx;
end
u={u1,u2,u3};
end