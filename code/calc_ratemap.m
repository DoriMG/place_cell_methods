function rms = calc_ratemap(binned_pos, n_bins, df_f) 

rms = nan(size(df_f,1), n_bins);
for i = 1:n_bins
    rms(:,i) = nanmean(df_f(:,binned_pos==i),2);
end