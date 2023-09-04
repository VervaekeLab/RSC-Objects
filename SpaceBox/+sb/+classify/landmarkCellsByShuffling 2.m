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
params.percentile = 99;
params.first_landmark_bins = 23:43;
params.second_landmark_bins = 63:83;
params.method = "mean"; % "mean": checks if the mean activity within the landmark bins window is above threshold; "individual": looks to see if there is a significant response in each individual bin for the landmark bins at both landmarks

% Print all variables if "help" is only input 
try
    if rmaps == "help"
        global verbose; verbose = true;
        sb.helper.displayParameters(params)
        verbose = false;
        return;
    end
end

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

    switch params.method
        case 'mean'

            % See if there is a response at first landmark position
            if nanmean(response(params.first_landmark_bins)) > nanmean(baseline_threshold(params.first_landmark_bins))
                response_first_position = 1;
            else
                response_first_position = 0;
            end

            % See if there is a response at second landmark position
            if nanmean(response(params.second_landmark_bins)) > nanmean(baseline_threshold(params.second_landmark_bins))
                response_second_position = 1;
            else
                response_second_position = 0;
            end

            significant_response = [response_first_position & response_second_position];

            if significant_response
               lcs = [lcs, c]; 
            end

        case 'individual'

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


end



