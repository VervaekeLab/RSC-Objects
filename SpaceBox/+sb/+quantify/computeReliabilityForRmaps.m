function reliability = computeReliabilityForRmaps(rmaps,ref_rmaps)

% rmaps: rastermap to compute reliability for each lap
% ref_rmaps: a reference rastermap that can be used to compute the baseline

% Parameters
reliability = [];   
std_required_for_significance = 3;
baseline_window_width = 10;
landmark_window = 30:60;

for c = 1:size(rmaps,3)

        % Find the n bins with the least activity
        roi_average = nanmean(rmaps(:,:,c));
        [~,sort_ind] = sort(roi_average);
        lowest_activity_inds = sort_ind(1:baseline_window_width);
        
        baseline = ref_rmaps(:,lowest_activity_inds,c);
        response = rmaps(:,landmark_window,c);
        
        baseline_mean = nanmean(baseline(:));
        baseline_std = nanstd(baseline(:));
                
        % For each lap, compare response to baseline
        lap_reliability = [];
        for l = 1:size(response,1)
            if nanmax(response(l,:)) > baseline_mean + (baseline_std * std_required_for_significance)
               lap_reliability(l) = 1;
            else
               lap_reliability(l) = 0;
            end
        end
        
        reliability(c) = sum(lap_reliability)/length(lap_reliability);
        
end

end