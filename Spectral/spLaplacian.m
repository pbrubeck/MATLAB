function del2 = spLaplacian(u, d)
% Calculates laplacian
del2=spPartialD(u, 1, 2);
for i=2:d
    del2=del2+spPartialD(u, i, 2);
end
end