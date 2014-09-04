

[horiz, min_range, max_range] = load_abenaki;
z_horiz = (horiz' - 2000)/4; % to match our 2000-3000ms subset

sz = size(z_horiz);
phase_map = NaN*ones(sz);
check_cube = zeros(sz);
for n = 1:sz(1)
	for m = 1:sz(2)
        z = ceil(z_horiz(n,m));
        if ~isnan(z)
            dd = phase_differences(n,m,z);
            while ( isnan(dd) && z < 251)
                z = z + 1;
                dd = phase_differences(n,m,z);
            end
            phase_map(n,m) = dd;
            check_cube(n,m,z) = 1;
        end
	end
end