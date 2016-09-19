function Dmf = sincdifft(f, M, h)

%  Dmf = sincdifft(f, M, h) computes the m-th derivative of the function
%  f(x) using the sinc differentiation process.   The function is 
%  assumed to be defined on the entire real line and the input values
%  correspond to samples of the function at N equispaced points
%  symmetric with respect to the origin.
%
%  Input:
%  f:    Vector of samples of f(x) at h*[-(N-1)/2:(N-1)/2]
%  M:    Number of derivatives required (integer).
%  h:    Step-size (real, positive).
%
%  Note:  0 < M < N-1.
%
%  Output:
%  Dmf:   m-th derivative of f
%
%  J.A.C. Weideman, S.C. Reddy 2000.   Corrected for complex data
%  by JACW, April 2003. 

f = f(:).';                 % Ensure f is a row vector
N = length(f);     
t = pi*[1:N-1];           

sigma = zeros(size(t));
                            % Compute first column & row of diff matrix
for l = 1:M;
    sigma = (-l*sigma + imag(exp(i*t)*i^l))./t;
      col = (pi/h)^l*[imag(i^(l+1))/(l+1) sigma];
end
      row = (-1)^M*col; row(1) = col(1);

% Imbed first row of Toeplitz matrix into bigger circulant matrix:
                            
rowbig = [row zeros(1,2^nextpow2(2*N)-2*N+1)  fliplr(col(2:N))];

% Multiply circulant matrix with data vector by using FFT:

    NN = length(rowbig);
     e = NN*ifft(rowbig);         % Eigenvalues of circulant matrix.
  fhat = fft([f zeros(1,NN-N)]);  % Take FFT of padded data vector,
   Dmf = ifft(e.*fhat);           % multiply the result by e-values and
   Dmf = Dmf(1:N).';              % take inverse FFT. 

if max(abs(imag(f))) == 0; Dmf = real(Dmf); end  % Real data in, real derivative out
