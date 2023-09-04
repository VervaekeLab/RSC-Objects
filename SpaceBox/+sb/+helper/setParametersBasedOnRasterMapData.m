function params = setParametersBasedOnRasterMapData(rmaps,params)
%% SET PARAMETERS BASED ON RASTERMAP DATA
% This is a helper function that takes a params struct and rastermap and
% sets the parameters default values based on the rastermap data.
%
% INPUT
%   rmaps: Rastermaps.
%   params: Struct containing key:value pairs. 
%   
% OUTPUT
%   params: Struct containing key:value pairs.
%
% Written by Andreas S Lande 2019


% Assume rastermap is of deconvolved signal if too small

if iscell(rmaps)
   rmaps = rmaps{1}; 
end

if max(rmaps(:)) < 0.5 
    params.color_map_minimum_max = 0.001;
end

end