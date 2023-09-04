function position_averaged_tuning = averagedTuningFromRasterMaps(rmaps,varargin)
%% averagedTuningFromRasterMaps
% Create a 1D averaged tuning map from a 2D rastermap. The function can
% also take a 3D input to generate a 2D rastermap where each row is the
% averaged tuning of rastermap n.
%
% INPUT
%   -- Required
%   rmaps: Matrix of dim: (Y x X x N) containing N rastermaps. rmaps can
%       also be a cell array with a 3D rmaps as each element, organized as
%       follow {rmaps1,rmaps2,...rmapsN}. If so, the maps are merged into a
%       single 2D rmap before given as output.
%
%   -- Optional
%   varargin: Alternating argument name and value, ex: ["keyword", value].
%       Arguments - See "Default parameters".
%
% OUTPUT
%   position_averaged_tuning: 1D array containing the averaged of the 2D
%       maps given as input.
%
% Written by Andreas S Lande 2019

%% Default parameters
params = struct();
params.smoothing = 1; % Smoothing applied to the 1D map after averaging.
params.smoothing_method = "Gaussian"; % Method of smoothing (Gaussian, movmean, etc .. see MATLAB function smoothdata()).
params.method = 'mean'; % Can be mean or median

% Update parameters
params = sb.helper.updateParameterStruct(params,varargin);

%% Create positionAveragedTuning
if iscell(rmaps) % Multiple maps to be merged
    all_rmaps = rmaps;
    counter = 1;
    
    % Init empty
    position_averaged_tuning = [];
    
    for x = 1:length(all_rmaps)
        
        % Get each rmaps from rmaps array
        rmaps = all_rmaps{x};
        
        % Average all maps within umaps
        for y = 1:size(rmaps,3)
            
            % Average rastermap
            switch params.method
                case 'mean'
                    position_averaged_tuning(counter,:) = nanmean(rmaps(:,:,y));
                case 'median'
                    position_averaged_tuning(counter,:) = nanmedian(rmaps(:,:,y));
            end
            
            counter = counter + 1; 
            
        end
        
    end
    
else
    
    % Average all rmaps
    for y = 1:size(rmaps,3)
        
        % Average rastermap
        switch params.method
            case 'mean'
                position_averaged_tuning(y,:) = nanmean(rmaps(:,:,y));
            case 'median'
                position_averaged_tuning(y,:) = nanmedian(rmaps(:,:,y));
        end
        
        % Apply smoothing
        if params.smoothing > 0
            position_averaged_tuning(y,:) = smoothdata(position_averaged_tuning(y,:),params.smoothing_method,params.smoothing);
        end
        
    end
    
end


end