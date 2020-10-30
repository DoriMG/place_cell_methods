function [data_norm]=find_delta_f(data2norm, frames2norm)
% Normalises data using the delta F over F method with F defined as a baseline period
% Adapted from code written by dr Kira Shaw 
%
% Inputs:
%   data2norm   - N x t array with fluorescence of N cells over time (t)
%   frames2norm - vector containing the frames to be included as the
%                 baseline


data_norm = nan(size(data2norm));
for i=1:size(data2norm,1)
    %subtract baseline, then divide by baseline
    data_norm(i,:)= (data2norm(i,:)-(nanmean(data2norm(i,frames2norm(1):frames2norm(2)))))...
        ./nanmean(data2norm(i,frames2norm(1):frames2norm(2)));
end

