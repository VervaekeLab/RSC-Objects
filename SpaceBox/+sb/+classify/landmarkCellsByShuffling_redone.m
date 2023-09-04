function lcs = landmarkCellsByShuffling(rmaps,varargin)
%% landmarkCellsByShuffling
%
% INPUT
%   -- Required
%
%
%   -- Optional
%   varargin: Alternating argument name and value, ex: ["keyword", value].
%       Arguments - See "Default parameters" section.
%
% OUTPUT
%
%


%% Default parameters
params = struct();
params.percentile = 95;
params.first_landmark_bins = 23:43;
params.second_landmark_bins = 63:83;

% Print all variables if "help" is only input 

% Update parameters
params = sb.helper.updateParameterStruct(params,varargin);
% Init empty var
lcs = [];

%% Shuffle data
for c = 1:size(rmaps,3)
    rmap = rmaps(:,:,c);
    response = nanmean(rmap);

    % Shuffle 1000 times
    shuffled_avg = zeros(1000,size(response,2));
    for i = 1:1000

        % Init zeroed matrix
        shuffled_rmap = zeros(size(rmap,1),size(rmap,2));

        % Shuffle for each lap
        for l = 1:size(rmap,1)
           shuffled_rmap(l,:) = circshift(rmap(l,:),randi(size(response,2))); 
        end

        shuffled_avg(i,:) = nanmean(shuffled_rmap);            

    end

    baseline_threshold = prctile(shuffled_avg,params.percentile);
    
    % First landmark
    response_first_position = [];
    for b = params.first_landmark_bins
        if response(b)>baseline_threshold(b)
            response_first_position = [response_first_position, b];
        end
    end
    
    % Second landmark
    response_second_position = [];
    for b = params.second_landmark_bins
        if response(b)>baseline_threshold(b)
            response_second_position = [response_second_position, b];
        end
    end
    
    
    if length(response_first_position)>0
        if length(response_second_position)>0
            lcs = [lcs,c];
        end
    end
    
    


end



