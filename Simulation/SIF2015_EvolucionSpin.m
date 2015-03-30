%%%%%%%%%%%%%%%%%%%% Taller SIF 2015 %%%%%%%%%%%%%%%%%%%%

% Author: Hector Sosa Martinez

% Rutina para calcular la evolucion del espin cuantico bajo la influencia de un campo magnetico 
% variante en el tiempo, para una particula con espin 1/2.

%%%%%%%%%%%%%%%%%%%% Taller SIF 2015 %%%%%%%%%%%%%%%%%%%%

clear all; close all; clc

%%%%%%% Parametros iniciales %%%%%%%

time = 2*pi;              % Total time for the time evolution
n = 100;                  % Intervals between initial and final time
t = linspace(0,time,n);   % Time vector

w_0 = 1;                  % Larmor frequency
dt = t(2)-t(1);           % Time step

% Diferent initial states
% state = [1;0];                    % Initial Spin State of the Atom (|+>_z)
% state = 1/sqrt(2)*[1;1];            % Initial Spin State of the Atom (|+>_x)
state = 1/sqrt(2)*[1;i];            % Initial Spin State of the Atom (|+>_y)


%%%%%%% Campo Magnetico %%%%%%%

A = 0.5;                    % Range = [0 1], Parameter A for the magnetic field
omega = 20*pi*0.01;       % Frecuency of the perpendicular magnetic field

u_x = A*cos(omega*t);         % x-component of the Magnetic field
u_y = A*sin(omega*t);         % y-component of the Magnetic field
u_z = sqrt(1-A^2)*ones(1,n);  % z-component of the Magnetic field


%%%%%%% Spin Evolution %%%%%%%

for ii = 1:n
    S_x(ii) = [conj(state(1,1,ii)),conj(state(2,1,ii))]*[0,1;1,0]*[state(1,1,ii);state(2,1,ii)];     % S_x
    S_y(ii) = [conj(state(1,1,ii)),conj(state(2,1,ii))]*[0,-1i;1i,0]*[state(1,1,ii);state(2,1,ii)];  % S_y
    S_z(ii) = [conj(state(1,1,ii)),conj(state(2,1,ii))]*[1,0;0,-1]*[state(1,1,ii);state(2,1,ii)];    % S_z
    S(ii) = sqrt(S_x(ii)^2+S_y(ii)^2+S_z(ii)^2);                                                     % S
    P_state(ii) = [conj(state(1,1,ii)),conj(state(2,1,ii))]*[state(1,1,ii);state(2,1,ii)];           % <s|s>
    
    H(:,:,ii) = w_0*[sqrt(1-A^2), A*exp(-1i*omega*t(ii));...     % Time dependent Hamiltonian
                     A*exp(1i*omega*t(ii)), -sqrt(1-A^2)];
    U(:,:,ii) = expm(-1i*H(:,:,ii)*dt);                          % "Infinitesimal" time evolution operator
    state(:,:,ii+1) = U(:,:,ii)*state(:,:,ii);                   % New Spin state after evolution
end

figure(1)  % Magnetic field over time
plot(t,u_x,'^b',t,u_y,'or',t,u_z,'+g')
legend('u_x','u_y','u_z')

figure(2)  % Expectation value of S over time
plot(t,S_x,'^b',t,S_y,'or',t,S_z,'+g')
legend('S_x','S_y','S_z')

figure(3)  % Validation quantities
plot(t,P_state,'b',t,S,'+r')
legend('<s|s>','|S|')


% %%%%%%% Visual representation %%%%%%%

% Creation of the sphere
    num = 50;
    eta = linspace(0,pi,num);
    fi = linspace(0,2*pi,num);

    x_1 = sin(pi/2).*cos(fi);
    y_1 = sin(pi/2).*sin(fi);
    z_1 = cos(pi/2)*ones(1,num);
    x_2 = sin(eta).*cos(0);
    y_2 = sin(eta).*sin(0);
    x_3 = sin(eta).*cos(pi/2);
    y_3 = sin(eta).*sin(pi/2);
    x_4 = sin(eta).*cos(pi);
    y_4 = sin(eta).*sin(pi);
    x_5 = sin(eta).*cos(3*pi/2);
    y_5 = sin(eta).*sin(3*pi/2);
    z_2 = cos(eta);

for ii = 1:n
    handle = figure(4);
    set(handle,'defaultlinelinewidth',1,'position',[100,100,512,512],'color','w','paperpositionmode','auto')
    subplot('position',[0,0,1,1]);
    hold on;
    plot3([1.2,-1.2],[0,0],[0,0],'k');
    plot3([0,0],[1.2,-1.2],[0,0],'k'); 
    plot3([0,0],[0,0],[1.2,-1.2],'k'); 
    plot3(x_1,y_1,z_1,'--k',x_2,y_2,z_2,'--k',x_3,y_3,z_2,'--k',x_4,y_4,z_2,'--k',...
          x_5,y_5,z_2,'--k')
    plot3([0,u_x(ii)],[0,u_y(ii)],[0,u_z(ii)],'b','LineWidth',3)
    plot3([0,S_x(ii)],[0,S_y(ii)],[0,S_z(ii)],'r','LineWidth',3)
    hold off;
    title('Spin Evolution')
    xlabel('x'); ylabel('y'); zlabel('z')
    legend('Magnetic Field in Blue','Spin Vector in Red')
    legend('boxoff')
    axis equal; axis off;
    view(136,28)

    h = getframe(gcf);
    Animacion(ii) = h;

    pause(.01)
    
    close(4);
end

% Creation of the movie

movie(Animacion,1,n);        
movie2avi(Animacion,'TimeEvol','fps',10,'quality',100);  