function [PC,p, rms, peak_value] = peak_method(df_f, location, run_frames, n_bins, frame_rate, n_reps, p_threshold)

% PEAK_METHOD  Find cells classified as place cells using the Peak method.
%   Inputs:
%   df_f:       Array containing the fluorescence of the ROIs in
%               delta F/F (ROIs x frames)
%   location:   Vector with the location (in cm) of the animal
%               at each frame (frames x 1)
%   run_frames: Boolean vector of frames in which animal is
%               running (frames x 1)
%   n_bins:     Number of location bins
%   frame_rate: Frame rate of the recording (Hz)
%   n_reps:     Number of shuffle repeats
%   p_threshold:Percentile threshold for cells to be included (optional)
%
%   Outputs:
%   PC:         Boolean vector of cells classified as place cells
%   p:          Probability of being a place cell for each cell (0-1)
%   rms:        Rate map for each cell
%   peak_value:     Peak value of each cell

if nargin <7
    p_threshold = 99;
end

n_rois = size(df_f,1);

% Bin the location trace
binned_pos = discretize(location(run_frames), n_bins);

% Calculate ratemap
rms = calc_ratemap(binned_pos, n_bins, df_f(:,run_frames));  
peak_value = squeeze(max(rms,[],2)); % Find peak value

% Assign array to store shuffle values
peak_shuffles = nan(size(df_f,1), n_reps);

%Perform the shuffles
for reps = 1:n_reps
    % Get random offset (> 5 seconds)
    shufMin = ceil(frame_rate*5);
    shufMax = size(df_f,2)-2*shufMin;
    randOffset = randi(shufMax)+shufMin;
    
    % Shuffle traces and calculate rms
    shuffTrace = circshift(df_f,randOffset,2);
    rmsshuff = calc_ratemap(binned_pos, n_bins, shuffTrace(:,run_frames));  
    
    % Calculate peaks of shuffled values
    peak_shuffles(:,reps) = max(rmsshuff,[],2);
end

% Find the threshold and place cell probability (p) for each cell
thresh = nan(n_rois,1);
p = nan(n_rois,1);
for roi = 1:n_rois
    thresh(roi) = prctile(peak_shuffles(roi, ~isnan(peak_shuffles(roi,:)) & ~isinf(peak_shuffles(roi,:))),p_threshold);
    p(roi) = invprctile(peak_shuffles(roi, ~isnan(peak_shuffles(roi,:)) & ~isinf(peak_shuffles(roi,:))), peak_value(roi))/100;
end

% Place cells are cells with a peak value higher than the threshold for
% that cell
PC =  peak_value>=thresh;

