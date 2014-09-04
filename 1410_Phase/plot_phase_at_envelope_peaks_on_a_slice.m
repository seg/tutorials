% Load Slice of SEGY
[data, segyHeader] = ReadSegy('data/penobscot_xl1155.sgy');

% Compute the complex attributes
[z, envelope, phase] = complex_attributes_on_section(data);

sz = size(data);
envelope_peaks = zeros(sz);
phase_at_envelope_peaks = zeros(sz);
idealised_phase = zeros(sz);
for n = 1:sz(2)
	% Detect peaks in envelope using non-maximal supression
	[envelope_peaks(:,n), phase_at_envelope_peaks(:,n)] = find_peaks(envelope(:,n), phase(:,n));
	idealised_phase(:,n) = get_idealised_phase(envelope_peaks(:,n), phase(:,n));
end

figure(4)
imagesc(fliplr(phase_at_envelope_peaks(500:750,:)))
colormap(hsv(256))
colorbar
set(gca, 'ytick', [1 50 100 150 200 250]);
set(gca, 'yticklabel', {'2000','2200','2400','2600','2800','3000'});
set(gca,'FontSize',axes_font_size)
title('Phase at Envelope Peaks', 'FontSize', title_font_size)


figure(5)
phase_error = abs(idealised_phase - phase_at_envelope_peaks);
phase_error(:,190) = 1.5;
imagesc(fliplr(phase_error(500:750, :)))
colormap(jet(256))
colorbar
set(gca, 'ytick', [1 50 100 150 200 250]);
set(gca, 'yticklabel', {'2000','2200','2400','2600','2800','3000'});
set(gca,'FontSize',axes_font_size)
title('Absolute Phase Error at Envelope Peaks', 'FontSize', title_font_size)
