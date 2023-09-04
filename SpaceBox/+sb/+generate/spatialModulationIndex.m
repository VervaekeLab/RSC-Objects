function [SMI, peak_index] = spatialModulationIndex(rmaps)
% Compute the Spatial modulation index for rastermaps of putative landmark
% cells.

method = 'ratio';

peak_index = [];

% For each cell
for c = 1:size(rmaps,3)
    
    
    switch method
        
        
        case 'ratio'
            % Compute the ratio between the two, so basically dividing
            % first and second landmark amplitude
            
            % Create position tuning map for the cell and normalize it
            avg_tuning = normalize(nanmean(rmaps(:,:,c)),'range',[0,1]);
            
            % Find the max position around 50 cm +/- 15 cm and at 110 +/- 15 cm.
            [~,first_peak] = max(avg_tuning(23:43));
            [~,second_peak] = max(avg_tuning(63:83));
            
            % Set correct position
            first_peak = first_peak + 23;
            second_peak = second_peak + 63;
            
            % Find the average response around the two peaks in a 10 cm window
            first_response = nanmean(avg_tuning(first_peak-3:first_peak+3));
            second_response = nanmean(avg_tuning(second_peak-3:second_peak+3));

            % Find biggest and smallest
            [biggest_response,biggest_ind] = max([first_response,second_response]);
            smallest_response = min([first_response,second_response]);
            
            SMI(c) = biggest_response/smallest_response;
            peak_index(c) = biggest_ind;
            
            if c == 62
               a = 1; 
            end
            
        case 'between_-1_and_1'
            
            % Create position tuning map for the cell and normalize it
            avg_tuning = normalize(nanmean(rmaps(:,:,c)),'range',[0,1]);
            
            % Find the max position around 50 cm +/- 15 cm and at 110 +/- 15 cm.
            [~,first_peak] = max(avg_tuning(23:43));
            [~,second_peak] = max(avg_tuning(63:83));
            
            % Set correct position
            first_peak = first_peak + 23;
            second_peak = second_peak + 63;
            
            % Find the average response around the two peaks in a 10 cm window
            first_response = nanmean(avg_tuning(first_peak-3:first_peak+3));
            second_response = nanmean(avg_tuning(second_peak-3:second_peak+3));
            
            % Compute SMI which is normalized between -1 and 1, where -1 means
            % the first landmark is full response, 0 means equal amplitude, and
            % 1 means second landmark has full response.
            SMI(c) = (second_response - first_response)/(second_response + first_response);
            
            % Visualize for debugging
            if false
                figure(1); clf;
                plot(avg_tuning,'Color','k','LineWidth',2);
                hold on;
                line([first_peak-3,first_peak+3],[0.8,0.8],'Color','r')
                line([second_peak-3,second_peak+3],[0.8,0.8],'Color','r')
                title(SMI(c))
            end
            
        case 'diamanti2020'
            % Compute SMI like Diamanti et al 2020.
            
            % Create average tuning map for even and odd laps
            avg_tuning_odd = normalize(nanmean(rmaps(1:2:end,:,c)),'range',[0,1]);
            avg_tuning_even = normalize(nanmean(rmaps(2:2:end,:,c)),'range',[0,1]);
            
            % Find the first and second peak position using odd laps map
            % Find the max position around 50 cm +/- 15 cm and at 110 +/- 15 cm.
            [~, first_peak] = max(avg_tuning_odd(23:43));
            [~, second_peak] = max(avg_tuning_odd(63:83));
            
            % Set correct position
            first_peak = first_peak + 23;
            second_peak = second_peak + 63;
            
            % Find the average response around the two peaks in a 10 cm window
            first_response_odd = nanmean(avg_tuning_odd(first_peak-3:first_peak+3));
            second_response_odd = nanmean(avg_tuning_odd(second_peak-3:second_peak+3));
            
            [~,preferred_response_odd] = max([first_response_odd,second_response_odd]);
            [~,not_preferred_response_odd] = min([first_response_odd,second_response_odd]);
            
            response_even = [];
            response_even(1) = nanmean(avg_tuning_even(first_peak-3:first_peak+3));
            response_even(2) = nanmean(avg_tuning_even(second_peak-3:second_peak+3));
            
            preferred_response_even = response_even(preferred_response_odd);
            not_preferred_response_even = response_even(not_preferred_response_odd);
            
            SMI(c) = (preferred_response_even - not_preferred_response_even) / (preferred_response_even + not_preferred_response_even);
            
    end
end



end