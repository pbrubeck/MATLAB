
%  The script file sineg.m solves the sine-Gordon equation
%  u_tt-u_xx-sin u on the real line using one of the following
%  differentiation matrices: (1) Hermite, (2) sinc, or (3) Fourier.
%  The solution is displayed as a mesh plot.

%  J.A.C. Weideman, S.C. Reddy 1998.  
%  Modified for MATLAB R13 by JACW in April 2003. 


method = input(' Which method: (1) Hermite, (2) sinc, (3) Fourier? ');

     N = input(' Order of differentiation matrix: N = ? ');

tfinal = input(' Final time: t = ? ');

if method == 1
b = input(' Scaling parameter for Hermite method: b = ? ');
[x,D] = herdif(N,2,b);               % Compute Hermite differentiation matrices
    D = D(:,:,2);                    % Extract second derivative
elseif method == 2
h = input(' Step-size for sinc method: h = ? ');
[x,D] = sincdif(N,2,h);              % Compute sinc differentiation matrices
    D = D(:,:,2);                    % Extract second derivative
elseif method == 3
L = input(' Half-period for Fourier method: L = ? ');
[x,D] = fourdif(N,2);                % Compute Fourier second derivative
    x = L*(x-pi)/pi;                 % Rescale [0, 2pi] to [-L,L]
    D = (pi/L)^2*D;
end 

     u0 = zeros(size(x));            % Compute initial conditions
     v0 = 2*sqrt(2)*sech(x/sqrt(2));   
     w0 = [u0; v0];

 tspan  = [0:tfinal/40:tfinal];
options = odeset('RelTol',1e-6,'AbsTol',1e-6);      % Options for ODE45
  [t,w] = ode45(@sgrhs, tspan, w0, options, D);     % Solve ODEs


u = w(:,1:N);                         % Extract u variable from solution array

mesh(x,t,u); view(30,30);             % Generate a mesh plot of u

