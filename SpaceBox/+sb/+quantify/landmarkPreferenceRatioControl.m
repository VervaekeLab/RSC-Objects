function ratio = landmarkPreferenceRatioControl(rmaps)
% Compute the landnmark preference but not by comparing the two
% landmark response amplitudes. Instead a ratio is computed on the same
% landmark but using even vs odd laps.

% Compute the ratio between the two, so basically dividing
% first and second landmark amplitude
i = 1;
for c = 1:size(rmaps,3)
    
    
    % Create position tuning map for the cell and normalize it
    avg_tuning_odd = normalize(nanmean(rmaps(1:2:end,:,c)),'range',[0,1]);
    avg_tuning_even = normalize(nanmean(rmaps(2:2:end,:,c)),'range',[0,1]);

    % Find the max position around 50 cm +/- 15 cm and at 110 +/- 15 cm.
    [~,first_peak] = max(avg_tuning_odd(23:43));
    [~,second_peak] = max(avg_tuning_odd(63:83));

    % Set correct position
    first_peak = first_peak + 23;
    second_peak = second_peak + 63;

    % Find the average response around the two peaks in a 10 cm window
    first_response_odd = nanmean(avg_tuning_odd(first_peak-3:first_peak+3));
    first_response_even = nanmean(avg_tuning_even(first_peak-3:first_peak+3));
    second_response_odd = nanmean(avg_tuning_odd(second_peak-3:second_peak+3));
    second_response_even = nanmean(avg_tuning_even(second_peak-3:second_peak+3));

    % Find biggest and smallest for looking at the first landmark
    [biggest_response,biggest_ind] = max([first_response_odd,first_response_even]);
    smallest_response = min([first_response_odd,first_response_even]);
    
    ratio(i) = ((biggest_response/smallest_response)-1) * 100;
    
    % If the first peak is the biggest, then set ratio to negative.
    if biggest_ind == 1
        ratio(i) = - ratio(i);
    end
    
    if isinf(ratio(i)) || (ratio(i)> 1000)
        ratio(i) =NaN;end
    i = i + 1;
%     
%     % Find biggest and smallest for looking at the second landmark
%     [biggest_response,biggest_ind] = max([second_response_odd,second_response_even]);
%     smallest_response = min([second_response_odd,second_response_even]);
%   
%     ratio(i) = ((biggest_response/smallest_response)-1) * 100;
%     
%     % If the first peak is the biggest, then set ratio to negative.
%     if biggest_ind == 1
%         ratio(i) = - ratio(i);
%     end
%     
%     i = i + 1;
%     
       
end


            
end