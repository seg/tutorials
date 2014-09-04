% Load a seimsic trace
x = sin(0:pi/4:4*pi);

% Compute the fft
X = fft(x, 2048);
Xphase = angle(X);

shifts = [0 30 90 180];
for n = 2:4
	x(n,:) = fftshifter(x(1,:), shifts(n)*pi/180);
	X(n,:) = fft(x(n,:), 2048);
	Xphase(n,:) = angle(X(n,:));
end

figure(2)
clf
for n = 1:4
	subplot(2,2,n)
	bar(real(X(n,:)))
	hold on
	set(bar(imag(X(n,:))), 'edgecolor', 'r')
	axis tight
	title(['Phase Shift ' num2str(shifts(n))])	
end