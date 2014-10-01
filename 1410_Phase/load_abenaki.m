function [H, min_range, max_range] = load_abenaki()
data = load('-ascii','data/abenaki_horizon.ascii');
min_range = min(data);
max_range = max(data);

%H = NaN*ones(max_range(1) - min_range(1), max_range(2) - min_range(2));
H = NaN*ones(601, 481);

sz = size(data);
for i = 1:sz(1)
	a = data(i,1) - min_range(1) + 1;
	b = data(i,2) - min_range(2) + 1;
	H(a, b) = data(i,3);
end