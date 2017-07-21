close all
clear

jdqr, pause

jdqr('ILU',2,30,8,struct('Disp',1)), pause

jdqr('ILU',2.1,30,8,struct('Disp',1)), pause

jdqr('Example1',2.1,30,8,struct('Disp',1),'ILU'), pause

jdqr('Example2',0.0,30,8,struct('Disp',1),'Example2'), pause

jdqr('Example2',0.0,30,8,struct('Disp',1),'ILU'), pause

jdqr('Example2',0.0,30,8,struct('Disp',1),'Example2','Example2'), pause


A=[1 2 3;2 4 5;1 9 8];  jdqr(A,struct('Disp',1)), pause
L=[1 0 0; 0 1 0; 0 0 1]; jdqr(A,struct('Disp',1),L), pause
M=[L,L]; jdqr(A,struct('Disp',1),M), pause

A=rand(20,20); jdqr(A,struct('Disp',1)), pause
