function [] = dumbbell(N, k)
a=1;
b=0.5;
L=pi;
vertex=[L+(a+1i*L)/2, (a+1i*L)/2, (a+1i*b)/2];
vertex=[vertex, fliplr(-conj(vertex))];
vertex=[vertex, fliplr(conj(vertex))];
corners=[1 6 7 12];
rectlap(vertex, corners, N, k);
end