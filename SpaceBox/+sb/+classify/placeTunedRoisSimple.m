function [place_tuned, stats] = placeTunedRoisSimple(rmaps,varargin)
%% placeTunedRoisSimple
% Using rastermaps for one or several rois, each rastermap is classified as
% either place tuned or not according to the algorithm below.
%
% INPUT
%   -- Required
%   rmaps: Matrix of dim: (Y x X x N) containing N roi rastermaps.
%
%   -- Optional
%   varargin: Alternating argument name and value, ex: ["keyword", value].
%       Arguments
%       + verbose (Default: 1): Int 0 or 1 to specify how much text output
%           should be generated during the execuiton of the function to update
%           user on the classification process.
%
% OUTPUT
%   place_tuned: Boolean array indicating whether the rastermap at each index
%       is place tuned or not. (dim: 1 x size(rmaps,3)) 
%   stats: Struct containing tuning statistics for each rastermap within
%       rmaps. (dim: 1 x size(rmaps,3)
%   
% Written by Andreas S Lande 2019

%% Default parameters
params.min_peak_dist = 40; % cm. Distance needed for two peaks to be considered as individual. This is X bins, not cm. Each bin is 1.5 cm typically.
params.min_peak_width = 2; % cm. Mao uses 15 cm, but for deconvolved signals of highly precise place cells that does not make sense, the aim is 2-3 cm here.
params.max_peak_width = 120; % cm.
params.max_n_peaks = 3; % Max number of peaks to be detected in the peak detection algorithm.
params.placefield_threshold_value = 0.30; 
params.reliability = 0.33; % 0.25 means 25% of trials / laps must have activity peak within detected field
params.activity_ratio_criteria = 3; % Arbitrary ratio criteria used also in Mao et al 2017
params.bin_size = 1.5; % Ex. bin size in cm
params.smoothing = 0; % Number of bins to be used for smoothing the rastermaps. If set to 0 smoothing is omitted.
params.smoothing_method = "gaussian"; % See smoothdata() MATLAB function for list of methods.

% Update parameters
params = sb.helper.updateParameterStruct(params,varargin);

%% Output parameters used if verbose > 0
if params.verbose
    fprintf('\n\n\n\n-- Classifying place tuned cells --\n');
    sb.helper.displayParameters(params);
end

%% Run classification
n_bins = size(rmaps,2);
stats = struct();
for roi = 1:size(rmaps,3)

    % Create positon averaged tuning map
    position_map = nanmean(rmaps(:,:,roi));

    % Smooth position_map
    if params.smoothing > 0
       position_map =  smoothdata(position_map,params.smoothing_method,params.smoothing);
    end
    
    stats(roi).position_map = position_map;
    stats(roi).rastermap = rmaps(:,:,roi);
    
    % Define place field threshold for given position map
    pf_threshold = min(position_map)+(params.placefield_threshold_value*(max(position_map)-min(position_map))); % Threshold is set as n % of the difference between the highest and lowest activity of the position activity map.
    
    % Add the beginning of the actvity map to the end to remove issues with not detecting peaks at edges
    position_map = [position_map position_map(1:30)];
    
    % Use the MATLAB Signal Processing Toolbox function findpeaks() to localize peaks within the position activity map
    [peaks,locs] = findpeaks(position_map,'MinPeakDistance',round(params.min_peak_dist/params.bin_size),'MinPeakHeight',pf_threshold,'WidthReference','halfheight','MinPeakWidth',round(params.min_peak_width/params.bin_size),'MaxPeakWidth',round(params.max_peak_width/params.bin_size),'NPeaks',params.max_n_peaks,'SortStr','descend');
    
    % Correct the bin number for the edge correction used above.
    locs(locs>n_bins) = locs(locs>n_bins)-n_bins;
    position_map = position_map(1:n_bins);
    locs = unique(locs);
    
    % Check to see that a single peak has not been mistaken as two at the
    % edge case
    try
        if length(peaks) == 2
            loc_shifted = abs(locs-n_bins);
            if abs(loc_shifted(1)-locs(2)) < 5
                peaks = peaks(1);
                locs = locs(1);
            end
        end
    end
    
    % If no peaks are found, add zeros to all stats fields.
    if length(locs) == 0
        
        % Set all values to zero
        stats(roi).intra_field_value = 0;
        stats(roi).field_width = 0;
        stats(roi).field_peak_ind = 0;
        stats(roi).field_start = 0;
        stats(roi).field_end = 0;
        stats(roi).n_peaks = 0;
        stats(roi).in_field_inds = 0;
        stats(roi).out_of_field_inds = 0;
        stats(roi).out_of_field_value = 0;
        stats(roi).activity_ratio = 0;
      
    else % Else run statistics analysis
        
        % For each peak detected
        for peak = 1:length(locs)
            peak_indx = locs(peak);
            
            % Find all neighbouring values that are above the threshold to obtain
            % place field width
            pf_start = peak_indx; % start index of place field
            pf_end = peak_indx; % end index of place field
            
            % Go from peak index backwards to locate beginning of place field
            look = 1;
            curr_pos = peak_indx;
            while look
                curr_pos = curr_pos - 1;
                if curr_pos < 1
                    curr_pos = n_bins;
                end
                
                if position_map(curr_pos) > pf_threshold
                    pf_start = curr_pos;
                else
                    look = 0;
                end
            end
            
            % Go forward from peak to find end index of place field
            look = 1;
            curr_pos = peak_indx;
            while look
                curr_pos = curr_pos + 1;
                if curr_pos > n_bins
                    curr_pos = 1;
                end
                
                if position_map(curr_pos) > pf_threshold
                    pf_end = curr_pos;
                else
                    look = 0;
                end
            end
            
            % Get all in-field and out-of-field indices
            in_field_inds = [];
            out_of_field_inds = [];
            if pf_end < pf_start
                in_field_inds =  [in_field_inds, 1:pf_end];
                in_field_inds = [in_field_inds, pf_start:length(position_map)];
                out_of_field_inds = [pf_end+1:pf_start-1];
            else
                in_field_inds = [in_field_inds, pf_start:pf_end];
                out_of_field_inds = [1:pf_start-1,pf_end+1:length(position_map)];
            end
            
            in_field_value = nanmean(position_map(in_field_inds));
            out_of_field_value = nanmean(position_map(out_of_field_inds));
            
            % Add place field information to stats
            stats(roi).intra_field_value(peak) = in_field_value;
            stats(roi).field_width(peak) = length(in_field_inds);
            stats(roi).field_peak_ind(peak) = locs(peak);
            stats(roi).field_start{peak} = pf_start;
            stats(roi).field_end{peak} = pf_end;
            stats(roi).n_peaks = length(locs);
            stats(roi).in_field_inds{peak} = in_field_inds;
            stats(roi).out_of_field_inds{peak} = out_of_field_inds;
            stats(roi).peaks = peaks;
            stats(roi).out_of_field_value(peak) = out_of_field_value;
            
        end
        
        % For each field, obtain out-of-field/in-field ratio
        total_out_of_field_inds = 1:n_bins; % Start with the assumption that all bins are out of field
        total_in_field_inds = [stats(roi).in_field_inds{:}];
        stats(roi).total_in_field_inds = total_in_field_inds;
        total_out_of_field_inds = total_out_of_field_inds(~ismember(total_out_of_field_inds,total_in_field_inds));
        stats(roi).total_out_of_field_inds = total_out_of_field_inds;
        for peak = 1:length(locs)
            stats(roi).out_of_field_value(peak) = nanmean(position_map(total_out_of_field_inds));
            stats(roi).activity_ratio(peak) = stats(roi).intra_field_value(peak)/stats(roi).out_of_field_value(peak);
        end  
    end
end
        
 
% Count number of times the peak activity for a lap is within the place
% field of a roi
for roi = 1:length(stats)
    
    n_laps_with_peak_in_field = 0;
    for lap = 1:size(rmaps,1)
        
        [max_val,ind] = max(stats(roi).rastermap(lap,:));
        
        if max_val == 0
            ind = 0;
        end
        
        % Check that the peak activity is not just noise
        if max_val > -inf
            
            % Is the peak activity of the lap within any of the place field peak positions?
            try
                in_field_inds = unique([stats(roi).in_field_inds{:}]);
                
                if ismember(ind,in_field_inds)
                    n_laps_with_peak_in_field = n_laps_with_peak_in_field + 1;
                end
                
            catch
                n_laps_with_peak_in_field = 0;
            end

        end
    end
    
    stats(roi).reliability = n_laps_with_peak_in_field/size(rmaps,1);
   
end

% Classify as place tuned if relibaly active and have a high activity ratio
for roi = 1:length(stats)
    
    % Start with assumption that roi is place tuned, set to zero if any
    % check below fails.
    stats(roi).isPlaceTuned = 1;

    % Check reliability
    if stats(roi).reliability < params.reliability
        stats(roi).isPlaceTuned = 0;
    end
    
    % Check activity ratio
    if stats(roi).activity_ratio < params.activity_ratio_criteria
        stats(roi).isPlaceTuned = 0;
    end
    
end

% Create a logical array with place tuned or not for all roi rastermaps
place_tuned = [stats(:).isPlaceTuned];

% Convert place_tuned into index instead of boolean array
place_tuned = find(place_tuned);

% Display number of place tuned rois if params.verbose == 1
if params.verbose
    fprintf('RESULT: %i out of %i (%.1f%%) rois place tuned (using placeTunedRoisSimple).\n',length(place_tuned),length([stats(:).isPlaceTuned]),length(place_tuned)/length([stats(:).isPlaceTuned])*100);
end
   
   
end