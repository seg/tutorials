function ideal_phase = get_idealised_phase(peaks, phase)
% flip to nearest 'zero' phase point
ideal_phase = zeros(size(phase));
ideal_phase(abs(phase) < pi/2) = 0;
ideal_phase(abs(phase - pi) < pi/2) = pi;
ideal_phase(abs(phase + pi) < pi/2) = -pi;

phase(peaks == 0) = NaN;
ideal_phase(peaks == 0) = NaN;