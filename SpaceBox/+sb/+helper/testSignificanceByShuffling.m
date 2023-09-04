function [significant_response, baseline_threshold] = testSignificanceByShuffling(rmap,bins,varargin)
%% testSignificanceByShuffling
%
% INPUT
%   -- Required
%       rmap: rastermap
%       bins: the x-axis bins of the rastermap that will be used to
%           determine significance
%
%   -- Optional
%   varargin: Alternating argument name and value, ex: ["keyword", value].
%       Arguments - See "Default parameters" section.
%
% OUTPUT
%   % significant_response : binary
%   % baseline_threshold: value at each x axis bin for the given percentile 

%% Default parameters
params = struct();
params.threshold = 99;
params.n_shuffles = 1000;
params.method = 'avg'; % avg: average of all rows in rmap. individual: significance test for each row in rmap.
params.plot = true;

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


%% Do shuffling

switch params.method
    
    case 'avg'
        response = nanmean(rmap);

        % Shuffle 1000 times
        shuffled_avg = zeros(params.n_shuffles,size(rmap,2));
        for i = 1:params.n_shuffles

            % Init zeroed matrix
            shuffled_rmap = zeros(size(rmap,1),size(rmap,2));

            % Shuffle for each lap
            for l = 1:size(rmap,1)
               shuffled_rmap(l,:) = circshift(rmap(l,:),randi(size(rmap,2))); 
            end

            shuffled_avg(i,:) = nanmean(shuffled_rmap);            

        end

        baseline_threshold = prctile(shuffled_avg,params.threshold);

        if nanmean(response(bins)) > nanmean(baseline_threshold(bins))
            significant_response = 1;
        else
            significant_response = 0;
        end
         
end

% Visualize
if params.plot
    
    
    figure(191991);
    clf;
    plot(response);
    hold on;
    plot(baseline_threshold,'--');
    
    
end









