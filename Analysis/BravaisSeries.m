function h = BravaisSeries(f, A, n)
% Returns the Fourier expansion of a Bravais-lattice-periodic-function.
[n1, n2, n3]=ndgrid(0:n-1);
B=A/n*[n1(1:end); n2(1:end); n3(1:end)];
h=reshape(f(B(1,:), B(2,:), B(3,:)), n, n, n);
h=fftn(h)/n^3;
end
