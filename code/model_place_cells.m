function [all_traces, all_loc] = model_place_cells(varargin)
% Creates a dataset containing model place cells and random control cells
%
% vargargin inputs:
%   bin_size        - size of the bins in cm (default: 2cm)
%   pf_width        - width of the place field in cm (default: 50 cm)
%   n_total         - total number of cells to be included in the dataset (default: 100)
%   perc_rand       - proportion of the dataset that are control cells (default: 0.8)
%   tot_n_tras      - number of traversals to include in the dataset (default: 50)
%   max_peak        - maximum fluorescence (in deltaF/F) at the peak of a place field (default: 1.2873)
%   var_sigma       - variance allowed in the location of the peak of the place field
%                     per traversal as a proportion of the place field width (default: 0)
%   rel             - reliability as the proportion of traversals each place cell fires in (default: 1)
%   env_length      - the length of the environment of the imported traversal in cm (default: 200 cm)
%   noise_lambda    - the lambda value used to create the Poisson noise (default: 235.1003 cm)

% outputs:
%   all_traces      - N x t array with the activity of N cells over time (t)
%   all_loc         - t x 1 vector containing the location at each timepoint (t) in
%                     cm

% Load pre-calculated parameters
all_tras = importdata(('allTras.mat'));

% Default values
bin_size = 2;
pf_width = 50;
n_total = 100;
perc_rand = 0.8;
tot_n_tras = 50;
n_placefields = 1;
max_peak = 1.2873;
var_sigma = 0;
rel  = 1;
env_length = 200;
noise_lambda = 235.1003;

% Parse input
p = inputParser;
addParameter(p,'bin_size',bin_size);
addParameter(p,'pf_width',pf_width);
addParameter(p,'n_total',n_total);
addParameter(p,'perc_rand',perc_rand);
addParameter(p,'tot_n_tras',tot_n_tras);
addParameter(p,'max_peak',max_peak);
addParameter(p,'var_sigma',var_sigma);
addParameter(p,'rel',rel);
addParameter(p,'n_placefields',n_placefields);
addParameter(p,'env_length',env_length);
addParameter(p,'noise_lambda',noise_lambda);
parse(p,varargin{:});

% put parsed values in parameters
bin_size = p.Results.bin_size;
pf_width = p.Results.pf_width;
n_total = p.Results.n_total;
perc_rand = p.Results.perc_rand;
tot_n_tras = p.Results.tot_n_tras;
max_peak = p.Results.max_peak;
var_sigma = p.Results.var_sigma;
rel = p.Results.rel;
n_placefields = p.Results.n_placefields;
env_length = p.Results.env_length;
noise_lambda = p.Results.noise_lambda;

% Calculate number of bins
n_bins = round(env_length/bin_size);

% Calculate number of place cells
perc_pcs = 1-perc_rand;
n_pcs = round(n_total*perc_pcs);

% Pick out random set of traversals and generate loco
rand_tras = datasample(1:length(all_tras),tot_n_tras);
all_loc = vertcat(all_tras{rand_tras});

% Clean loco
all_loc = all_loc-min(all_loc);
all_loc = all_loc/max(all_loc);
raw_loc = all_loc*(env_length-1)+1;
all_loc = round(raw_loc);

% Calculate the average location of the place cells
offsets = linspace((pf_width/2)*n_placefields, env_length-(pf_width/2)*n_placefields, n_pcs+1)-(env_length/2);
pf_centres = [(pf_width/2):pf_width:pf_width*n_placefields];
pf_centres = pf_centres-nanmean(pf_centres)+(n_bins/2);

% Create place cells
traces = []; % Save eventual traces
x = [1:bin_size:env_length]; % Location vector
for i = 1:n_pcs
    % Determine in what traversals the cell should fire (depends on
    % reliability (rel)
    travs_on = datasample(1:tot_n_tras,rel*tot_n_tras,'Replace',false);
    fire_trav = zeros(tot_n_tras,1);
    fire_trav(travs_on) = 1;

    temp_trace = [];
    % Model fluorescence per traversal
    for j = 1:length(rand_tras)
        trace = zeros(size(all_tras{rand_tras(j)}));
        if fire_trav(j)==1
            % Determine the center of firing in this traversals (depends on
            % variability (var_sigma)
            center = normrnd(pf_centres+offsets(i),var_sigma*pf_width);
            y = zeros(size(x));

            % Model fluorescence for each place field using a normal probability
            % density function
            for c = 1:length(center)
                y = y+normpdf(x,center(c),pf_width/4)';
            end
            y = y/max(y)*max_peak; % scale the peak
            temp_loc = all_tras{rand_tras(j)};


            for k = 1:length(trace)
                trace(k) = y(round(temp_loc(k)/bin_size)); %TODO: Dori wtf?
            end
        end

        temp_trace = [temp_trace; trace]; % append trace for each traversal
    end
    traces = [traces, temp_trace]; % append all traces
end

% Add random traces
rand_traces = poissrnd(noise_lambda, n_total, length(all_loc));
rand_traces = find_delta_f(rand_traces,1:100);
all_traces = [traces'; zeros(n_total*perc_rand, length(all_loc)) ];
all_traces = all_traces+rand_traces;
