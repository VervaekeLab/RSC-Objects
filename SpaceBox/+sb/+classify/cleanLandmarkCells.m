function [lcs, tcs] = cleanLandmarkCells(rmaps,rmaps2, varargin)
%% cleanLandmarkCells
% Some landmark cells have very clear landmark responses but they are not
% detected by the traditional algorithm because either of the bumps are too
% weak. Yet, looking at their rastermaps they have clearly responses at
% both presentations of the landmark. This function tries to classify these
% cells, by looking for cells that have activity close to the landmark but
% almost a averaged zero response everywhere else.
% To exclude the obvious false positives, all rastermaps are sorted by
% their peak position and only cells that have a max response around the
% landmarks are included.
%
% INPUT
%   -- Required
%      rmaps: Rastermap of cells
%
%   -- Optional
%   varargin: Alternating argument name and value, ex: ["keyword", value].
%       Arguments - See "Default parameters" section.
%
% OUTPUT
%   lcs: index of all clean lcs.
%   tcs: index of all clean persisting landmark cells.


%% Default parameters
params = struct();

% Print all variables if "help" is only input
try
    if inp == "help"
        global verbose; verbose = true;
        sb.helper.displayParameters(params)
        verbose = false;
        return;
    end
end

% Update parameters
params = sb.helper.updateParameterStruct(params,varargin);


%% Sort the average responses of these cells
[sorted_avg,sort_ind] = sb.sort.rasterMapByPeak(sb.generate.averagedTuningFromRasterMaps(rmaps));

%% Classify
lcs = [];
tcs = [];

for c = 1:size(rmaps,3)
    
    rmap = rmaps(:,:,c);
    threshold = nanmean(rmap(:,[1:30,90:105]),1);
    min_peak_height = prctile(threshold,99);
    
    [~,peaks] = findpeaks(nanmean(rmap),'NPeaks',20,'MinPeakDistance',10,'MinPeakHeight',min_peak_height,'SortStr','descend');
     
    lcs(c) = 0;
    
    if length(peaks) >= 2
        peaks = peaks(1:2);    
        
        if ismember(peaks(1),32:50)
            if ismember(peaks(2),72:90)
                
                % Check that there is a minimum reliability
                realibility = 0;
                
                
                
                if reliability > 0.25
                    lcs(c) = 1;
                end
            end
        end
    end
    
    
    % Check if the cell is also a persisting landmark cell
    if(size(rmaps,3) == size(rmaps2,3))
        rmap2 = rmaps2(:,:,c);
        [~,peaks2] = findpeaks(sb.generate.averagedTuningFromRasterMaps(rmap2),'NPeaks',2,'MinPeakHeight',nanmean(rmap2(:)));

        tcs(c) = 0;
        if length(peaks2) == 2

            if ismember(peaks2(1),32:38)
                if ismember(peaks2(2),72:78)
                    tcs(c) = 1;
                end
            end
        end

    end
    
    
end


lcs = find(lcs);
tcs = find(tcs);




