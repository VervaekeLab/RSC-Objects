function active_rois = activeRoisAPT(dff,varargin)
%% activeRoisAPT
% Find active ROIs by applying activity peak thresholding (APT). Every signal
% in dff will be checked to see that max of the signal is above the
% params.activity_peak_threhold value.
%
% INPUT
%   -- Required
%   dff: Neuropil subtracted dFF signals for all rois of interest. Shape of
%       input is M x N, where M is number of rois, and N is number of samples.
%  
%   -- Optional
%   varargin: Alternating argument name and value, ex: ["keyword", value].
%       Arguments
%       + activity_peak_threshold (Default: 0.8): The minimum peak activity
%           of the dff signal.
%
% OUTPUT
%   active_rois: Boolean array indicating in each given index whether the
%       roi in dff is considered "active" or not.
%
% EXAMPLE
%   1) Finding all ROIs with minimum dFF peak above 1.5 is done as follows:
%       active_rois = sb.classify.activeRoisAPT(dff,'activity_peak_threshold',1.5);
%
% Written by Andreas S Lande 2019

%% Default parameters
params = struct();
params.activity_peak_threshold = 0.3;

% Update parameters
params = sb.helper.updateParameterStruct(params,varargin);

%% Display info
if params.verbose
    fprintf('\n-- Classyfing active cells --\n');
end

%% Run classification
% Set all rois to default not active
active_rois = zeros(1,size(dff,1));

% For each roi, set active if max activity is greater than
% threshold value
for x = 1:size(dff,1)
    if max(dff(x,:))> params.activity_peak_threshold
        active_rois(x) = 1;
    end
end

% Display info
if params.verbose > 0 
    fprintf('RESULT: %i of %i (%.1f%%) rois active (using activeRoisAPT)\n',sum(active_rois),length(active_rois),(sum(active_rois)/length(active_rois))*100);
end

% Convert active rois to indexing instead of boolean array
active_rois = find(active_rois);
end