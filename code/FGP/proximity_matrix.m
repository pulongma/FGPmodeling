function W = proximity_matrix(locs, d1)




n = size(locs, 1);
nzmx = 20*n; 
W = spalloc(n, n, nzmx);
for i=1:n
	dtemp = pdist2(locs(i, :), locs);
	indi = find(dtemp<d1);
	W(indi, i) = 1;
	W(i, i) = 0;
end


