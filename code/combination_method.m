function [PC,p]= combination_method(df_f, location,velocity, varargin)

% COMBINATION_METHOD  Find cells classified as place cells using the
% Combination method. Adapted from Dombeck (2010).
%   Inputs:
%   df_f:       Array containing the fluorescence of the ROIs in
%               delta F/F (ROIs x frames)
%   location:   Vector with the location (in cm) of the animal
%               at each frame (frames x 1)
%   velocity:   Vector with speed of animal (in cm/s) (frames x 1)
%   varargin:
%       speedThresh (2) - The minimum speed (cm/sec) for inclusion
%       runLength (40) - The minimum continuous run length (cm)
%       segmentsThresh (20) - The minimum number of runs needed
%       trackLength (200) - The length of the track (cm)
%       nbins (100) - The number of bins to divide the track into
%       smth (3) - The width of the boxcar smoothing kernel
%       baselineP (0.1333) - The percentile threshold for the basline
%           florescence
%       inFieldThresh (0.25) - The fraciontal distance between the baseline
%           and peak florescence to be used as a threshold for "in field"\
%       minPeakHeight (0.1) - At least one point in an accepted field must
%           exceed this value
%       inOutRatioThresh (4) - The in-to-out of field intensity ratio must
%           exceed this to be accepted
%       minFieldWidth (20) - The minimum width of an accepted place field
%       maxFieldWidth (150) - The maximum width of an accepted place field
%       activeRunsThresh (0.2) - The minimum ratio of runs in which the
%       cell should be active
%       p_threshold (95) - Percentile threshold for cells to be included
%
%   Outputs:
%   PC:         Boolean vector of cells classified as place cells
%   p:          Probability of being a place cell for each cell (0-1)
%
%   Adapted from:
% 1.	Sheffield, M. E. J. & Dombeck, D. A. Calcium transient prevalence
%       across the dendritic arbour predicts place field properties. Nature
%       517, 200–204 (2015).
% 2.	Sheffield, M. E. J., Adoff, M. D. & Dombeck, D. A. Increased
%       Prevalence of Calcium Transients across the Dendritic Arbor during
%       Place Field Formation. Neuron 96, 490–504.e5 (2017).
% 3.	Dombeck, D. A., Harvey, C. D., Tian, L., Looger, L. L. & Tank, D.
%       W. Functional imaging of hippocampal place cells at cellular
%       resolution during virtual navigation. Nat. Neurosci. 13, (2010).

ip = inputParser;
ip.addParamValue('speedThresh',2);
ip.addParamValue('runLength',40);
ip.addParamValue('segmentsThresh',20);
ip.addParamValue('nbins',100);
ip.addParamValue('smth',3);
ip.addParamValue('baselineP',20/150);
ip.addParamValue('inFieldThresh',0.25);
ip.addParamValue('minPeakHeight',0.1);
ip.addParamValue('p_threshold',95);
ip.addParamValue('inOutRatioThresh',4);
ip.addParamValue('activeRunsThresh',0.2);
ip.addParamValue('minFieldWidth',20);
ip.addParamValue('bootstrap',1e3);
ip.addParamValue('trackLength',200);
ip.addParamValue('maxFieldWidth',120);
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end

% Check there are enough running epochs
epochs = longRunningEpochs(location,velocity,'speedThresh',speedThresh,'runLength',runLength);
if size(epochs,1)<segmentsThresh
    warning(sprintf('Not enough long-running segments: %i of %i',size(epochs,1),segmentsThresh));
    
    % Return empty output
    if size(df_f,2)>1
        flds = cell(size(df_f,2),1);
    else
        flds = [];
    end
    PC = zeros(1,size(df_f,2));
    p = inf(1,size(df_f,2));
    return
end


% Build occupancy map & loc list
% Make a logical vector of frames to include
vid_goods = zeros(size(location));
for i = 1:size(epochs,1)
    vid_goods(epochs(i,1):epochs(i,2)) = 1;
end

goods = find(vid_goods);

% Make a histogram of bin dwell times (oc) and the bin visited at each
% time
binEdges = linspace(0,trackLength,nbins+1);
[oc,loc] = histc(location(goods),binEdges);


% Build spatial intensity maps
sim = arrayfun(@(i)nanmean(df_f(goods(loc==i),:),1),1:nbins,'UniformOutput',false);
sim = cat(1,sim{:});

% pim = arrayfun(@(i)nansum(Fc3(goods(loc==i),:)>0,1)/sum(loc==i),1:nbins,'UniformOutput',false);
% pim = cat(1,pim{:});

% Boxcar smoothing
sim = conv2(full(sim),padarray(ones(smth,1)/smth,[1 1],0,'both'),'same');

% Find the baseline for each cell
baseline = sort(sim);
baseline=mean(baseline(1:round(nbins*baselineP),:));

% Find the potential fields
potentialFields = logical(zeros(size(sim)));
for i = 1:size(sim,1)
    potentialFields(i,:) = sim(i,:)>(max(sim)-baseline)*inFieldThresh;
end

% Find their bounds
ons = arrayfun(@(i)threshEpochs(potentialFields(:,i)),1:size(df_f,2),'UniformOutput',false);
offs = cellfun(@(ons)ons(:,2),ons,'UniformOutput',false);
ons = cellfun(@(ons)ons(:,1),ons,'UniformOutput',false);


% 1.) Field width thresholds
goods = cellfun(@(on,off)(off-on)*mean(diff(binEdges)),ons,offs,'UniformOutput',false);
goods = cellfun(@(x)x>minFieldWidth&x<maxFieldWidth,goods,'UniformOutput',false);

ons = cellfun(@(x,y)x(y),ons,goods,'UniformOutput',false);
offs = cellfun(@(x,y)x(y),offs,goods,'UniformOutput',false);

% 2.) The field must have one value of at least minPeakHeight mean DF/F
goods = arrayfun(@(i)arrayfun(@(j)max(sim(ons{i}(j):offs{i}(j),i))>minPeakHeight,(1:numel(ons{i}))'),1:numel(ons),'UniformOutput',false);
peakHeights = arrayfun(@(i)arrayfun(@(j)max(sim(ons{i}(j):offs{i}(j),i)),(1:numel(ons{i}))'),1:numel(ons),'UniformOutput',false);

ons = cellfun(@(x,y)x(y),ons,goods,'UniformOutput',false);
offs = cellfun(@(x,y)x(y),offs,goods,'UniformOutput',false);

% 3.) The mean in field DF/F value must be inOutRatioThresh times the mean out of field DF/F value
temp = (1:nbins);
% Make a logical vector for each cell labeling the potential in-field
% regions
for f = 1:size(ons,2)
    if isempty(ons{f})
        inField{f} = false(size(sim,1),1);
    else
        inField{f} = any(temp>ons{f}&temp<offs{f},1);
    end
end


goods = arrayfun(@(i)arrayfun(@(j)nanmean(sim(ons{i}(j):offs{i}(j),i))/...in field
    nanmean(sim(~inField{i},i))...out of field
    >inOutRatioThresh,(1:numel(ons{i}))'),1:numel(ons),'UniformOutput',false);
inOutField = cellfun(@max,arrayfun(@(i)arrayfun(@(j)nanmean(sim(ons{i}(j):offs{i}(j),i))/...in field
    nanmean(sim(~inField{i},i))...out of field
    ,(1:numel(ons{i}))'),1:numel(ons),'UniformOutput',false),'UniformOutput',false);

ons = cellfun(@(x,y)x(y),ons,goods,'UniformOutput',false);
offs = cellfun(@(x,y)x(y),offs,goods,'UniformOutput',false);

% 4.) The cell must be active during the whole time the animal is in the
goods = {};
for i = 1:numel(ons)
    if numel(ons{i})
        for j=1:numel(ons{i})
            inField = location>=binEdges(ons{i}(j))&location<=binEdges(offs{i}(j));
            
            entrances = find(diff(inField)==1)+1;
            exits = find(diff(inField)==-1);
            if inField(1)
                entrances = [1;entrances];
            end
            if inField(end)
                exits = [exits;numel(inField)];
            end
            laps = [entrances exits];
            laps = laps(diff(laps,[],2)>1,:);
            
            exits = (1:numel(location))';
            entrances = NaN(size(laps,1),1);
            for k=1:numel(entrances)
                inField = exits>laps(k,1)&exits<laps(k,2)&vid_goods;
                if any(inField)
                    entrances(k) = any(df_f(inField,i)>0);
                end
            end
            goods{i}(j) = nanmean(entrances)>activeRunsThresh;
        end
    else
        goods{i} = [];
    end
end
ons = cellfun(@(x,y)x(y),ons,goods,'UniformOutput',false);
offs = cellfun(@(x,y)x(y),offs,goods,'UniformOutput',false);

% Turn to inidiviual field bounds & initialize output
% flds = cellfun(@(x,y)ifThenElse(isempty(x),@()[],[x(:) y(:)]),ons,offs,'UniformOutput',false);
for f = 1:size(ons,2)
    if isempty(ons{f})
        flds{f} = [];
    else
        flds{f} = [ons{f}(:) offs{f}(:)];
    end
end
isplace=cellfun(@(x)size(x,1),flds)>0;
p=zeros(size(isplace));
p(~isplace(:))=1;


if bootstrap>0% Bootstrapping necessary
    
    % Remove bootstrap from inputs so we don't recursivley bootstrap
    if any(cellfun(@(x)isequal(x,'bootstrap'),varargin))
        i = find(cellfun(@(x)isequal(x,'bootstrap'),varargin));
        varargin = varargin([1:i-1 i+2:end]);
    end
    
    nchunks = 5.8343e+09; %% Depends on size of computer memory
    nchunks=min(floor(nchunks/8/20/numel(df_f)),bootstrap);
    
    for j=find(isplace(:))'% Only bootstrap the candidate place cells
        %     for j = 1:length(isplace)
        
        try
            % Identify the discrete transients
            transients = regionprops(full(df_f(:,j)~=0),'PixelIdxList','Area');
            transients = cellfun(@(i)df_f(i),{transients.PixelIdxList},'UniformOutput',false);
            
            beta = sum(df_f(:,j)==0)/(numel(transients)+2);% Average gap length
            for CHUNK = 1:ceil(bootstrap/nchunks)
                
                bootstrap_ = min(bootstrap,CHUNK*nchunks)-(CHUNK-1)*nchunks;
                
                [~,order] = sort(rand(numel(transients),bootstrap_));% Randomly shuffle the order
                
                % Initialize gaps
                gaps = zeros(numel(transients)+2,bootstrap_);
                goods = false(1,bootstrap_);
                reps = 0;
                while any(~goods)
                    % Generate gaps from the exponential distribution
                    gaps(:,~goods) = exprnd(beta,numel(transients)+2,sum(~goods));
                    
                    % Force gaps to sum to the same length
                    gaps = round(gaps.*(sum(df_f(:,j)==0)./sum(gaps)));
                    
                    % See if we need to remake any
                    goods = sum(gaps)==sum(df_f(:,j)==0);
                    reps = reps+1;
                    if reps > 20
                        break
                    end
                end
                
                if reps <= 20
                    arfun = arrayfun(@(n)zeros(n,1),gaps(2:end-1,:),'UniformOutput',false);
                    if size(arfun,1) ~= size(transients(order),1)
                        arfun = arfun';
                    end
                    
                    Fc3_ = cat(3,transients(order),arfun);
                    Fc3_ = permute(Fc3_,[3 1 2]);
                    Fc3_ = reshape(Fc3_,[2*numel(transients) bootstrap_]);
                    Fc3_ = cat(1,arrayfun(@(n)zeros(n,1),gaps(1,:),'UniformOutput',false),Fc3_,arrayfun(@(n)zeros(n,1),gaps(end,:),'UniformOutput',false));
                    Fc3_ = reshape(cat(1,Fc3_{:}),size(df_f,1),bootstrap_);
                    
                    [~,temp,~]=placeTest(Fc3_,location,velocity,'bootstrap',0,varargin{:});% Turn off bootstrapping and test
                    p(j)=p(j)+sum(temp);
                else
                    p(j) = 1;
                end
            end
            
            p(j)=p(j)/bootstrap;% Make p value
        catch
            p(j) = 1;
        end
    end
    PC=p<((100-p_threshold)/100);
end


