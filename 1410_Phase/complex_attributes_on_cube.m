function [z, envelope, phase] = complex_attributes_on_cube(cube)
% function [z, envelope, phase] = complex_attributes_on_cube(section)
%
% expects data in (x,y,z) order
%

sz = size(cube);
hilbert = zeros(sz); 	% initialise an empty slice
for x = 1:sz(1)
	for y = 1:sz(2)			% for every trace
		% phase shift by 90 degrees so that the hilbert signal lags
		% the original trace. Hence -90, and resulting in generally
		% increasing phase values down the trace
		hilbert(x,y,:) = fftshifter(cube(x,y,:), -90); 
	end
end
z = cube + 1i*hilbert; % form a complex version of the section
envelope = sqrt(real(z).^2 + imag(z).^2);
phase = -1*angle(z);