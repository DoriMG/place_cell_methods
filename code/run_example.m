% Variables to set
n_bins = 100; % number of bins
env_length = 200; % length of environment in cm
bin_size = env_length/n_bins; % bin size in cm
frame_rate = 7.51; % frame rate in Hz
n_reps = 100; % number of shuffles to perform

% Create model populations
[df_f, all_loc] = model_place_cells('tot_n_tras', 100, 'perc_rand', 0.4);

% Preprocess loc trace and extract velocity and traversals
[loc_cm, vel, run_frames, traversals, trav_frame_ind] = preprocess_location(all_loc, env_length, 1, frame_rate); 


%% Run methods
% Peak method
[PC,p, rms, maxInt] = peak_method(df_f, all_loc, run_frames, n_bins, frame_rate, n_reps, 99);

%Information method
[PC,p, rms, si_value] = information_method(df_f, all_loc, run_frames, n_bins, frame_rate, n_reps, 95);

%Stability method
[PC,p, rms, stab_value] = stability_method(df_f, all_loc, run_frames, n_bins, trav_frame_ind, n_reps, 95);

%Combination method
[PC,p] = combination_method(df_f', all_loc,vel,'p_threshold',95,'segmentsThresh', 4);
