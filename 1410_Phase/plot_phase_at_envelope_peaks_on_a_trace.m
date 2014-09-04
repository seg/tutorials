% Load a seismic trace
% This is trace was exported from OpenDTect as a Simple File
[s, t] = load_simple_trace('data/penobscot_trace_il1190_xl1155.trace');
trace_length = numel(t);

% Compute the Analytic (Complex) Trace
z = s + i*ffthilbert(s);
envelope = abs(z);
phase = -1*angle(z);

% Detect peaks in envelope using non-maximal supression
[envelope_peaks, phase_at_envelope_peaks] = find_peaks(envelope, phase);
idealised_phase = get_idealised_phase(envelope_peaks, phase);
phase_differences = abs(idealised_phase - phase_at_envelope_peaks);

% plot
plot_range = 500:750;
tshort = t(plot_range);
h = figure(1);
clf
subplot(2,1,1)
plot(tshort, real(z(plot_range)),'k')
hold on
grid on
plot(tshort, envelope(plot_range),'r')
envelope_peaks(envelope_peaks == 0) = NaN;
stem(tshort, envelope_peaks(plot_range),'r');
title('Envelope', 'FontSize', title_font_size);
xlabel('Time, t(ms)', 'FontSize', axes_font_size);
set(legend('Original','Envelope','Envelope Peaks','location','northeast'),'FontSize', axes_font_size);
plot(tshort, -envelope(plot_range),'r')
stem(tshort, -envelope_peaks(plot_range),'r');
axis tight

subplot(2,1,2)
phase_at_envelope_peaks(envelope_peaks == NaN) = NaN;
stem(tshort, idealised_phase(plot_range), 'k');
hold on
grid on
stem(tshort, phase_at_envelope_peaks(plot_range),'r');
set(gca,'XLim',[tshort(1) tshort(end)]);

ax = findall(h,'type','axes');

title('Phase', 'FontSize', title_font_size)
xlabel('Time, t(ms)', 'FontSize', axes_font_size)
ylabel('Phase (Radians)', 'FontSize', axes_font_size)
set(legend('Actual Phase', 'Idealised Phase'),'FontSize', axes_font_size);