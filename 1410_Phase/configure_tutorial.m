current_dir = pwd;
addpath(current_dir);
addpath([current_dir '/SegyMAT']);
addpath([current_dir '/SegyMAT/GUI']);
disp('PATH updated');

% other config
pkg rebuild -noauto oct2mat
if (~isunix() && ~ismac())
	graphics_toolkit('fltk');
end
warning('off', 'Octave:possible-matlab-short-circuit-operator');

screensize = get(0,'screensize');
big_plot(1) = 0.1*screensize(3);
big_plot(2) = 0.1*screensize(4);
big_plot(3) = 0.9*screensize(3);
big_plot(4) = 0.9*screensize(4);

axes_font_size = 12;
title_font_size = 20;
