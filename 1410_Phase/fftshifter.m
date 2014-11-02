function y = fftshifter(x, phase_shift_in_radians)
% Create an array to apply the coefficient to 
% positive and negative frequencies appropriately
N = length(x);
M = ceil((N+1)/2);
R0 = exp(i*phase_shift_in_radians);
R = ones(size(x));
R(1:M) = R0;
R(M+1:N) = conj(R0);

% Apply the phase shift in the frequency domain
Xshifted = R.*fft(x);

% Recover the shifted time domain signal
y = real(ifft(Xshifted));