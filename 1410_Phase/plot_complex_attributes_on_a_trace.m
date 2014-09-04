% Load a seismic trace
% This is trace was exported from OpenDTect as a Simple File
[s, t] = load_simple_trace('data/penobscot_trace_il1190_xl1155.trace');

% Compute the Analytic (Complex) Trace
z = s + i*ffthilbert(s);
envelope = abs(z);
phase = angle(z);

% plot the complex trace
plot_range = 500:750; % select 1s window for plotting
tshort = t(plot_range);
h = figure(1)
clf
subplot(2,1,1)
plot(tshort, real(z(plot_range)),'k');
hold on;
plot(tshort, imag(z(plot_range)),'k-.');
plot(tshort, envelope(plot_range),'r');
plot(tshort, -envelope(plot_range),'r');
axis tight
set(gca,'FontSize', axes_font_size)
title('Trace Amplitudes', 'FontSize', title_font_size)
xlabel('Time, t(ms)', 'FontSize', axes_font_size)
set(legend('Original','Hilbert','Envelope','location','northeast'),'FontSize', axes_font_size)

subplot(2,1,2)
plot(tshort, phase(plot_range),'k');
axis tight
set(gca,'FontSize',axes_font_size)
title('Instaneous Phase')
xlabel('Time, t(ms)', 'FontSize', axes_font_size)
ylabel('Angle (radians)', 'FontSize', axes_font_size)
title('Original', 'FontSize', title_font_size);