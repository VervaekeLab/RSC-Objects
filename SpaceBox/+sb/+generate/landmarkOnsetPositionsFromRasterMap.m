function [x,y] = landmarkOnsetPositionsFromRasterMap(landmark_rmap,varargin)
%% landmarkOnsetPositionsFromRasterMap
% Find the onsets in x and y of the rastermap of where the landmarks are
% on.
%
% INPUT
%   -- Required
%       landmark_rmap: Rastermap showing the speed of the motor to be used.
%
%   -- Optional
%   varargin: Alternating argument name and value, ex: ["keyword", value].
%       Arguments - See "Default parameters" section.
%
% OUTPUT
%       onsets: Position of x and y for each landmark onset.
%
% Written by Andreas S Lande 2020


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


%% Find onset positions

% Create binary map
landmark_rmap(landmark_rmap>0) = 1;


y = [];
x = [];
for lap = 1:size(landmark_rmap,1)
    
    [~,p] = findpeaks(landmark_rmap(lap,:));
    
    y = [y,ones(1,length(p)) * lap];
    x = [x, p];
    
    
end
