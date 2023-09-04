function significant_response = detectSignificantResponsesByShuffling(rmap,bins,varargin)
%% detectSignificantResponsesByShuffling
% INPUT
%   rmap: (matrix) Rastermap to detect significance
%   bins: (array) bins that are used when comparing actual response to
%   shuffled
%
% OUTPUT
%   significant_response: true/false.


%% Default parameters
params = struct();
params.percentile = 95;

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


%% Shuffle data
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

if nanmean(response(bins)) > nanmean(baseline_threshold(bins))
    significant_response = 1;
else
    significant_response = 0;
end


end