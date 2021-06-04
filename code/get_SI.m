function all_si = get_SI(rm)

all_si = nan(size(rm,1),1);
for c = 1:size(rm,1)
    cell_temp = rm(c,:);
    cell_temp = cell_temp-(min(cell_temp));
    
    fi = cell_temp;
    f = nanmean(cell_temp);
    % Already normalized for occupancy
    
    mi = 0;
    for i = 1:size(rm,2)
        mi_temp = fi(i)*log2(fi(i)./f);
        if ~isnan(mi_temp)
            mi = mi+mi_temp;
        end
    end
    all_si(c) = mi;
end