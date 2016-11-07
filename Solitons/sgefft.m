function [Q,x] = sgefft(N, dt)
% Sine-Gordon Equation, fft linear propagator
% u_tt - u_xx = 0
x=sqrt(2*pi/N)*(-N/2:N/2-1)';
k=1i*sqrt(2*pi/N)*([0:N/2-1, -N/2:-1]');
q1=exp(dt*k);
q2=dt*sinc(1i*dt*k);
q3=exp(-dt*k);

function [u,v]=sgefftprop(u,v)    
    uhat=fft(u);
    vhat=fft(v);
    u=ifft(q1.*uhat+q2.*vhat,'symmetric');
    v=ifft(q3.*vhat,'symmetric');
end
Q=@sgefftprop;
end

