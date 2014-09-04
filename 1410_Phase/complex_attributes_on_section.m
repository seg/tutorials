function [z, envelope, phase] = complex_attributes_on_section(section)
% function [z, envelope, phase] = complex_attributes_on_section(section)

sz = size(section)
hilbert = zeros(sz); 	% initialise an empty slice
for n = 1:sz(2)			% for every trace
	% phase shift by 90 degrees so that the hilbert signal lags
	% the original trace. Hence -90, and resulting in generally
	% increasing phase values down the trace
	hilbert(:,n) = fftshifter(section(:,n), -90); 
end
z = section + 1i*hilbert; % form a complex version of the section
envelope = sqrt(real(z).^2 + imag(z).^2);
phase = -1*angle(z);


