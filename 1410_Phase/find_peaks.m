function [y, ya] = find_peaks(x, a)
% Detect peaks in envelope using non-maximal supression
window_length = 3;
trace_length = numel(x);
halflength = ceil((window_length - 1) / 2);
y = zeros(size(x));
ya = NaN*ones(size(x));
for n = int32(halflength+1:trace_length-halflength)
	range = n-halflength:n+halflength;
	local_max = max(x(range));
	if (x(n) >= local_max)
		y(n) = x(n);
		ya(n) = a(n);
	end
end