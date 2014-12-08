% Load Slice of SEGY
[data, segyHeader] = ReadSegy('data/penobscot_xl1155.sgy');

[z, envelope, phase] = complex_attributes_on_section(data);

plot_range = 500:750; % select 1s window for plotting

figure(3)
colormap(gray(256))
z(:,190) = 6500 + 1i*8000;
imagesc(fliplr(real(z(plot_range,:)))) % the original seismic section
set(gca, 'ytick', [1 50 100 150 200 250]);
set(gca, 'yticklabel', {'2000','2200','2400','2600','2800','3000'});
set(gca,'FontSize', axes_font_size)
title('Original', 'FontSize', title_font_size);
colorbar

figure(4)
imagesc(fliplr(imag(z(plot_range,:)))) % the (trace-wise) Hilbert transform section
set(gca, 'ytick', [1 50 100 150 200 250]);
set(gca, 'yticklabel', {'2000','2200','2400','2600','2800','3000'});
set(gca,'FontSize',axes_font_size)
title('Hilbert', 'FontSize', title_font_size)
colorbar

figure(5)
envelope(:,190) = 10000;
imagesc(fliplr(envelope(plot_range,:))) % the Envelope (Instantaneous Amplitude) = abs(z)
set(gca, 'ytick', [1 50 100 150 200 250]);
set(gca, 'yticklabel', {'2000','2200','2400','2600','2800','3000'});
set(gca,'FontSize',axes_font_size)
title('Envelope', 'FontSize', title_font_size)
colorbar

figure(6)
phase(:,190) = 0;
imagesc(fliplr(phase(plot_range,:))) % the Instantaneous Phase = angle(z)
set(gca, 'ytick', [1 50 100 150 200 250]);
set(gca, 'yticklabel', {'2000','2200','2400','2600','2800','3000'});
set(gca,'FontSize',axes_font_size)
title('Phase', 'FontSize', title_font_size)
colorbar