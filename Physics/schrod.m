function [] = schrod(n,m)
% Quantum harmonic oscillator time evolution in one spatial dimension
% A neat example on Hermite spectral methods

% Differentiantion matrix, nodes and quadrature weights
[D, x, w]=hermD(n+1);
dt=6.6/32^2;


% Quantum mechanical operators
H=-D*D+diag(x.^2);         % Halmintonian
A=D+diag(x);               % Annihilation operator
Q=expm(-1i*dt*H);          % Unitary time evolution operator
F=hermFT(x,w);             % Fourier Transform operator

% Eigenfunctions of A
[Psi, lambda]=eig(A); 
lambda=diag(lambda);
[lambda,order]=sort(lambda);
Psi=Psi(:, order(m));


% Normalization and interpolation to a finer, equispaced grid
w=w.*exp(x'.^2);
Psi=bsxfun(@rdivide, Psi, sqrt(w*(abs(Psi).^2)));
xx=linspace(x(1), x(end), 2048); xx=xx(:);
uu=interp1(x,Psi,xx,'spline');
vv=interp1(x,F*Psi,xx,'spline');


% Plot corresponding state
figure(1); clf;
subplot(3,1,1); 
h1=plot(xx, [real(uu), imag(uu)]);
ax=axis(gca); ax(1:2)=x([1 end]); ax(3:4)=[-1 1]*max(abs(ax(3:4))); axis(ax);
title('\Psi(x,t)');

subplot(3,1,2);
h2=plot(xx, abs(uu).^2); 
ax=axis(gca); ax(1:2)=x([1 end]); axis(ax);
title('|\Psi(x,t)|^2');

subplot(3,1,3);
h3=plot(xx, [real(vv), imag(vv)]); 
ax=axis(gca); ax(1:2)=x([1 end]); ax(3:4)=[-1 1]*max(abs(ax(3:4))); axis(ax);
title('\Psi(k,t)');

% Time evolution
for i=1:50000
    Psi=Q*Psi;
    uu=interp1(x,Psi,xx,'spline');
    vv=interp1(x,F*Psi,xx,'spline');
    
    set(h1, {'YData'}, {real(uu); imag(uu)});
    set(h2, 'YData', abs(uu).^2);
    set(h3, {'YData'}, {real(vv); imag(vv)});
    drawnow;
end
end

