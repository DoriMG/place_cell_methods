function [all_loc, vel, run_frames, traversals, run_frame_ind ] =preprocess_location(all_loc, env_length, convert, frame_rate); 

run_frame_thresh = 2;

%Convert location into cm
if convert
    all_loc = all_loc-min(all_loc);
    norm_loc = all_loc/max(all_loc);
    all_loc = norm_loc * env_length;
end

% Calculate velocity in cm
vel = [0;diff(all_loc)]*frame_rate;

% Find the traversals in the data (based on sudden change in location)
traversals = find_traversals(all_loc, 0.5*max(all_loc));
run_frames = vel>=2;

% Find the start and stop frames of individual traversals
startFrames = [1,find(diff(traversals)==1)];
stopFrames = [startFrames(2:end)];
stopFrames = [stopFrames, length(all_loc)];
run_frame_ind = [startFrames;stopFrames]';