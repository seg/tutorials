function [s, t, trace_header] = load_simple_trace(filename)
% This function will oad a trace that was exported from 
% OpenDTect as a Simple File with header info
s = load('data/penobscot_trace_il1190_xl1155.trace');
s = s(3:end); % extract data part and reverse to correct order
t = 0:4:6000;
trace_header = s(1:2);