NLSEThe following generates the paper's graphs.nlse(512, 0, 5, 20, 'bright', 'cheb')nlse(512, 0, 5, 20, 'dark2', 'cheb')nlse(512, -3, 3, 3, 'peregrine', 'cheb')The general form is nlse(N, t1, t2, L, init, method)Modify N to change number of collocation points. t1,t2 is time window, L is the space window plot in the interval [-L,L]; method corresponds to 'cheb' 'fft or 'herm'.

SGE

The following generates the paper's graphs.

sge(256, 0, 8, 4, 'breather', 'cheb')

sge(256, 0, 20, 4, 'kink', 'cheb')

The general form is 

sge(N, t1, t2, L, init, method)Modify N to change number of collocation points. t1,t2 is time window, L is the space window plot in the interval [-L,L]; method corresponds to 'cheb' 'fft or 'herm'.

