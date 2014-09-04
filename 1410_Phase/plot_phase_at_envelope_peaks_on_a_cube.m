% Load the penobscot3d cube
% load penobscot3D_loaded.mat

cube = reshape(data, 251, 481, 601);
cube = permute(cube, [2 3 1]);

% Compute the complex attributes
[z, envelope, phase] = complex_attributes_on_cube(cube);

sz = size(cube);
envelope_peaks = zeros(sz);
phase_at_envelope_peaks = zeros(sz);
idealised_phase = zeros(sz);
for n = 1:sz(1)
	for m = 1:sz(2)
		% Detect peaks in envelope using non-maximal supression
		[envelope_peaks(n,m,:), phase_at_envelope_peaks(n,m,:)] = find_peaks(envelope(n,m,:), phase(n,m,:));
		idealised_phase(n,m,:) = get_idealised_phase(envelope_peaks(n,m,:), phase(n,m,:));
	end
end

phase_differences = abs(idealised_phase - phase_at_envelope_peaks);

pd_cube = permute(phase_differences, [3 1 2]);
pd_cube = reshape(pd_cube, 251, 481 * 601);

% WriteSegyStructure('phase_diff_cube_fp9.sgy',SegyHeader,SegyTraceHeaders, int32(pd_cube));

% figure(4)
% imagesc(fliplr(phase_at_envelope_peaks(500:750,:)))
% colormap(hsv(256))
% colorbar
% set(gca, 'ytick', [1 50 100 150 200 250]);
% set(gca, 'yticklabel', {'2000','2200','2400','2600','2800','3000'});
% set(gca,'FontSize',axes_font_size)
% title('Phase at Envelope Peaks', 'FontSize', title_font_size)


% figure(5)
% phase_error = abs(idealised_phase - phase_at_envelope_peaks);
% phase_error(:,190) = 1.5;
% imagesc(fliplr(phase_error(500:750, :)))
% colormap(jet(256))
% colorbar
% set(gca, 'ytick', [1 50 100 150 200 250]);
% set(gca, 'yticklabel', {'2000','2200','2400','2600','2800','3000'});
% set(gca,'FontSize',axes_font_size)
% title('Absolute Phase Error at Envelope Peaks', 'FontSize', title_font_size)
