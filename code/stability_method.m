function [PC,p, rms, stability]= stability_method(df_f, location,run_frames, n_bins, run_frame_ind, n_reps,p_threshold)

% STABILITY_METHOD  Find cells classified as place cells using the Peak
% method. Adapted from O'Leary (2019).
%   Inputs:
%   df_f:       Array containing the fluorescence of the ROIs in
%               delta F/F (ROIs x frames)
%   location:   Vector with the location (in cm) of the animal
%               at each frame (frames x 1)
%   run_frames: Boolean vector of frames in which animal is
%               running (frames x 1)
%   n_bins:     Number of location bins
%   run_frame_ind:
%               Start and stop frames of each traversal (traversals x 2)
%   n_reps:     Number of shuffle repeats
%   p_threshold:Percentile threshold for cells to be included (optional)
%
%   Outputs:
%   PC:         Boolean vector of cells classified as place cells
%   p:          Probability of being a place cell for each cell (0-1)
%   rms:        Rate map for each cell
%   stability:  Peak value of each cell

if nargin <7
    p_threshold = 95;
end

n_rois = size(df_f,1);
n_runs = size(run_frame_ind,1);
run_frames = find(run_frames);

% Find first and second half of session
half_runs=round(n_runs/2);
runs_first_half = rangeConcat(run_frame_ind(1:half_runs,:));
runs_first_half = intersect(runs_first_half,run_frames);

runs_second_half = rangeConcat(run_frame_ind(half_runs+1:n_runs,:));
runs_second_half = intersect(runs_second_half,run_frames);

%bin the location data
location_first_half = discretize(location(runs_first_half), n_bins);
location_second_half = discretize(location(runs_second_half), n_bins);


%Calculate ratemaps for first and second half
rms_first_half = calc_ratemap(location_first_half, n_bins, df_f(:,runs_first_half));
rms_second_half = calc_ratemap(location_second_half, n_bins, df_f(:,runs_second_half));

%Calculate within cell and shuffled correlations
stability = nan(n_rois,1);
within_sesh_spat_corrs_roi_shuff = nan(n_rois,n_reps);
for r=1:n_rois
    %Calculate within cell correlation
    exc = isnan(rms_first_half(r,:))|isnan(rms_second_half(r,:)); %Exclude nans
    stability(r)= corr(rms_first_half(r,~exc)',rms_second_half(r,~exc)');
    
    %Perform shuffle
    for i = 1:n_reps
        rand_roi = randperm(n_rois,1); %Select random ROI
        while rand_roi == r % Make sure ROI is not compared to itself
            rand_roi = randperm(n_rois,1);
        end
        exc = isnan(rms_first_half(r,:))|isnan(rms_second_half(rand_roi,:)); %Exclude nans
        within_sesh_spat_corrs_roi_shuff(r,i) = corr(rms_first_half(r,~exc)',rms_second_half(rand_roi,~exc)');
    end
end

shuff_comp = reshape(within_sesh_spat_corrs_roi_shuff, 1,[]); % Flatten array

% Find the threshold and place cell probability (p) for each cell
thresh=prctile(shuff_comp,p_threshold);
p = nan(n_rois,1);
for n = 1:n_rois
    p(n) = invprctile(shuff_comp, stability(n))/100;
end

% Place cells are cells with a stability higher than the threshold for
% that cell
PC = stability>=thresh;
rms = cat(3,rms_first_half,rms_second_half); % save out ratemaps