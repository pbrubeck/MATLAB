function del2 = fftLap(u, d)
% Calculates laplacian
del2=fftD(u, 1, 2);
for i=2:d
    del2=del2+fftD(u, i, 2);
end
end