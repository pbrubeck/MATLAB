function [] = blackhole( N )
%BHM  Black Hole Merger
R1=10;
R2=5;
L=100;
th=0;
ww=bipolarMap(R1,R2,L,N);
figure(1);
h=mesh(real(ww),imag(ww),zeros(size(ww)));
c=4*(R1+R2+L); 
axis equal; view(2);
colormap([0,0,0]);

for i=1:100
    ww=exp(1i*th)*bipolarMap(R1,R2,L,N);
    L=L*0.98;
    th=th+pi/20;
    set(h, 'XData', real(ww));
    set(h, 'YData', imag(ww));
    xlim([-c,c]); ylim([-c,c]);
    drawnow; pause(0.1);
end
end