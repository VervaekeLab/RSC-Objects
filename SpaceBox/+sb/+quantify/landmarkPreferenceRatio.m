function ratio = landmarkPreferenceRatio(rmaps,varargin)
% Compute the ratio of preference for each cell between the two landmark
% positions. If the first landmark is twice as big as the second, it will
% be given a score of -100%, indicating a 100% increase difference towards 
% where - means it is the first landmark. 


if ~isempty(varargin)
    params = varargin{1};
else
    params.first_landmark_bins = 23:43;
    params.second_landmark_bins = 63:83;
end

% Compute the ratio between the two, so basically dividing
% first and second landmark amplitude
ratio = [];


for c = 1:size(rmaps,3)
    
    % Create position tuning map for the cell and normalize it
    avg_tuning = normalize(nanmean(rmaps(:,:,c)),'range',[0,1]);

    % Find the max position around 50 cm +/- 15 cm and at 110 +/- 15 cm.
    [~,first_peak] = max(avg_tuning(params.first_landmark_bins));
    [~,second_peak] = max(avg_tuning(params.second_landmark_bins));

    % Set correct position
    first_peak = first_peak + params.first_landmark_bins(1);
    second_peak = second_peak + params.second_landmark_bins(1);

    % Find the average response around the two peaks in a 10 cm window
    first_response = nanmean(avg_tuning(first_peak-3:first_peak+3));
    second_response = nanmean(avg_tuning(second_peak-3:second_peak+3));

    % Find biggest and smallest
    [biggest_response,biggest_ind] = max([first_response,second_response]);
    smallest_response = min([first_response,second_response]);

    ratio(c) = ((biggest_response/smallest_response)-1) * 100;
    
    % If the first peak is the biggest, then set ratio to negative.
    if biggest_ind == 1
        ratio(c) = - ratio(c);
    end
    
    % this removes strong outliers, these are falsifying actual
    % distribution and are a result of the shuffling method to detect
    % landmark cells
    if isinf(ratio(c))||ratio(c)> 1000 
        ratio(c) =NaN;end
       
end

end