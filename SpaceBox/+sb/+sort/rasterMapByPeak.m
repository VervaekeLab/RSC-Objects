function [sorted_map,sorted_peak_index] = rasterMapByPeak(rmap,varargin)
%% rasterMapByPeak
% Take as input a 2D raster map and sort each line in y dim by the max peak
% of each line.
%
% INPUT
%   -- Required
%   rmap: 2D raster map of M x N, where the sorting will be applied to
%       dimension M. 
%
%   -- Optional
%   varargin: Alternating argument name and value, ex: ["keyword", value].
%       Arguments - See "Default parameters".
%   
% OUTPUT
%   sorted_map: A map of similar size as rmap, where each row is now sorted
%       so that the first row is the row that has its peak earliest in the
%       column dimension of rmap.
%   sorted_peak_index: The index number of each element in rmap after
%       sorting. I.e. if the 5th element is first when sorted and the 7th
%       element is second after sorting, this variable would be [5,7, ...]
%
% Written by Andreas S Lande 2019

%% Default parameters
params = struct();
params.normalize = true; % Normalize each row in rmap to itself
params.normalize_type = 'range';
params.sort_order = 'peak'; % If peak it finds all peaks, else indices can be given as input and the output will be ordered accordingly.

% Update parameters
params = sb.helper.updateParameterStruct(params,varargin);

%% Normalize each row of raster map
if params.normalize
    switch params.normalize_type
        case 'range'
            for x = 1:size(rmap,1)
                 rmap(x,:) = normalize(rmap(x,:),'range');
            end        
            
        case 'zscore'
            rmap = normalize(rmap,2,'zscore');
    end
end

%% Sort rmap
if ischar(params.sort_order)
    if params.sort_order == 'peak'
        % Find the peak for each row
        [~,peak_ind] = nanmax(rmap');
        [~,sorted_peak_index] = sort(peak_ind);

        % Init sorted_map with all zeros
        sorted_map = zeros(size(rmap));

        % Sort the peak map matrix accordingly to the max peaks
        for p = 1:length(sorted_peak_index)
            sorted_map(p,:) = rmap(sorted_peak_index(p),:);
        end
    end
    
else
    
    
    % Sort the peak map matrix accordingly to the max peaks
    for p = 1:length(params.sort_order)
        sorted_map(p,:) = rmap(params.sort_order(p),:);
    end
    
    
    
    
end





