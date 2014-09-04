[horizon, min_range, max_range] = load_abenaki();
data = ReadSegyFast('data/penobscot3d_2400_2800ms.sgy');

X = reshape(data, 101, 481, 601);
X = permute(X, [2 3 1]);
clear data

sz = size(X);

horizon_mask = zeros(sz);
z = round((horizon - min_range(3)) / 4)';
z(z < 1) = 1;
z(z > sz(3)) = sz(3);

valid_z = ~isnan(z);
for n = 1:sz(1)
	for m = 1:sz(2)	
		if (valid_z(n,m))
			X(n,m, z(n,m)) = 1.5e9;
		end
	end
end

z_mean = ( max(max(z)) - min(min(z)) ) / 2;
z_shift = z - z_mean; 

flat = zeros(size(X));
for n = 1:sz(1)
    for m = 1:sz(2)
        if (valid_z(n,m))
            first = max([1 1 + z(n,m)]);
            last  = min([sz(3) sz(3) + z(n,m) - 10]);        
            flat(n,m,1:numel(first:last)) = X(n,m,first:last);
        end        
    end
end

imagesc(fliplr(squeeze(X(115,:,:))'))
imagesc(fliplr(squeeze(flat(115,:,:))'))




