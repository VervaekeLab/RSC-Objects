function [lcs,tcs,lcs_stats,tcs_stats] = landmarkAndTraceCells(rmaps,rmaps2,varargin)
% LANDMARK CELL METHOD 2
% Specify what landmark is removed, and thus check only for firing at this
% position still for rmaps2

% Use deconv.
% In this method I do the following:
% (1) Classify cells as LC based on having two peaks around the landmark
% positions of the motors. This is done on motor on laps only.
% Cells are only included if they have 2 peaks;
% the distance between these peaks has to be between 49 and 71 cm; and the
% two peaks has to be located at positon between 14 and 92 cm on wheel.
% (2) Each putative LC cell is then tested according to Mao criteria for
% place tuning.
% (3) Each cell then tested to see if the landmark that has been omitted
% still has a response in the averaged tuning map of the cell at +/- 20 cm
% from the landmark position.
% (4) The response at the omitted landmark position when motor is off has
% to be present on at least 20% of laps.

% Trace cells (tcs) are a subset of landmark cells (lcs) which is a subset
% of place tuned cells (pcs).

params.visualize = false;
params.detection_method = 'new';

% Update parameters
params = sb.helper.updateParameterStruct(params,varargin);

%% Find landmark cells based on rmaps
% Init
lcs_stats = struct();
lcs = [];
c = 1;

for r = 1:size(rmaps,3)
    
    % Create lap-averaged map of roi
    avg_map = nanmean(rmaps(:,:,r));
    
    % Smooth the averaged map
    avg_map_smoothed = smoothdata(avg_map,'Gaussian',20);
    
    % Set parameters for find peaks
    % Set a minimum threshold for a peak to be included
    peak_threshold = min(avg_map_smoothed)+(0.3*(max(avg_map_smoothed)-min(avg_map_smoothed))); % Threshold is set as 20 % of the difference between the highest and lowest activity of the position activity map.
    
    max_peak_width = 50;
    min_peak_dist = 25;
    n_peaks = 4;
    min_peak_width = 8;
    
    [~,peaks] = findpeaks(avg_map_smoothed,'NPeaks',n_peaks,'MinPeakHeight',peak_threshold,'MinPeakDistance',min_peak_dist,'SortStr','descend','MaxPeakWidth',max_peak_width,'MinPeakWidth',min_peak_width);
    
    % Use only two biggest peaks
    try
        peaks = peaks(1:2);
    end
    
    % Visualize if asked
    
    if params.visualize
        figure(1); clf;
        subplot(1,4,1:3);
        findpeaks(avg_map_smoothed,'NPeaks',n_peaks,'MinPeakHeight',peak_threshold,'MinPeakDistance',min_peak_dist,'SortStr','descend','MaxPeakWidth',max_peak_width,'MinPeakWidth',min_peak_width);
        hold on;
        
        line([33,33],[-1000,1000]);
        line([73,73],[-1000,1000]);
        line([93,93],[-1000,1000]);
        line([0,105],[peak_threshold,peak_threshold],'Color','r','LineStyle','-','LineWidth',2);
        %line([0,105],[nanmean(avg_map_smoothed),nanmean(avg_map_smoothed)],'Color','k','LineStyle','-');
        plot(avg_map,'Color','k','LineStyle','--','LineWidth',2);
        
        xtix = get(gca, 'XTick');
        set(gca, 'XTick',xtix, 'XTickLabel',xtix*1.5)
        xticks([0,33,73,105])
        xticklabels({'0','50','110','157'})
        xlabel('Position (cm)');
        set(gca,'FontSize',22);
        
        
        subplot(1,4,4);
        
        %imagesc([rmaps(:,:,r);rmaps2(:,:,r)]);
        imagesc(rmaps(:,:,r));
        title('NO');
        1;
        xtix = get(gca, 'XTick');
        set(gca, 'XTick',xtix, 'XTickLabel',xtix*1.5)
        xticks([0,33,73,105])
        xticklabels({'0','50','110','157'})
        xlabel('Position (cm)');
        set(gca,'FontSize',22);
        
    end
    
    % Only include cells with two peaks
    if length(peaks) == 2
        
        % Calculate peak to peak distance
        peak_distance = (max(peaks)-min(peaks))*1.5;
        
        % Check for ideal peak to peak distance
        if (peak_distance > 49) & (peak_distance < 71)
            
            % Check that peaks are not in beginning of end of track
            if (min(peaks) > 14) & (max(peaks) < 92)
                
                % Check that there is a low number bins where the avg_map
                % is above zero, to prune it down to only "sharp cells"
                % ...
                if sum(avg_map>peak_threshold)/length(avg_map) < 0.5
                    
                    lcs_stats(c).peaks = peaks;
                    lcs_stats(c).separation_distance = peak_distance;
                    c = c+1;
                    lcs = [lcs,r];
                    
                    if params.visualize
                        title('YES');
                    end
                end
            end
        end
    end
end


if length(lcs)>0
    
    % Find out which of the putative lcs are place tuned
    [pcs,pcs_stats] = sb.classify.placeTunedRoisSimple(rmaps(:,:,lcs),'reliability',0.25); % 15:92
    lcs = lcs(find(pcs));
    not_lcs = ones(1,size(rmaps,3));
    not_lcs(lcs) = 0;
    not_lcs = find(not_lcs);
    
    %% Find trace cells based on rmaps2
    % Init
    tcs_stats = struct();
    tcs = [];
    c = 1;
    
    % What bin is the landmark that has been manipulated
    manipulated_landmark_position = 33;
    
    switch params.detection_method
        
        case 'old'
            
            for r = 1:size(rmaps,3)
                
                % Create lap-averaged map of roi
                avg_map = nanmean(rmaps2(:,:,r));
                
                % Smooth the averaged map
                avg_map_smoothed = smoothdata(avg_map,'Gaussian',20);
                
                % Set a minimum threshold for a peak to be included
                peak_threshold = min(avg_map_smoothed)+(0.33*(max(avg_map_smoothed)-min(avg_map_smoothed))); % Threshold is set as 20 % of the difference between the highest and lowest activity of the position activity map.
                
                [~,peaks] = findpeaks(avg_map_smoothed,'NPeaks',4,'MinPeakHeight',peak_threshold,'MinPeakDistance',25,'SortStr','descend');
                
                % Use only two biggest peaks
                try
                    peaks = peaks(1:2);
                end
               
                
                % Visualize if asked
                %params.visualize = 1;
                if params.visualize
                    figure(1); clf;
                    subplot(1,4,1:3);
                    findpeaks(avg_map_smoothed,'NPeaks',4,'MinPeakHeight',peak_threshold,'MinPeakDistance',25,'SortStr','descend');
                    hold on;
                    
                    line([33,33],[-1000,1000]);
                    line([73,73],[-1000,1000]);
                    line([93,93],[-1000,1000]);
                    line([0,105],[peak_threshold,peak_threshold],'Color','r','LineStyle','-','LineWidth',2);
                    %line([0,105],[mean(avg_map_smoothed)+1*std(avg_map_smoothed),mean(avg_map_smoothed)+1*std(avg_map_smoothed)],'Color','r','LineStyle','--');
                    xtix = get(gca, 'XTick');
                    set(gca, 'XTick',xtix, 'XTickLabel',xtix*1.5)
                    xticks([0,33,73,105])
                    xticklabels({'0','50','110','157'})
                    xlabel('Position (cm)');
                    set(gca,'FontSize',22);
                    
                    subplot(1,4,4);
                    %imagesc([rmaps(:,:,r);rmaps2(:,:,r)]);
                    imagesc(rmaps2(:,:,r));
                    xtix = get(gca, 'XTick');
                    set(gca, 'XTick',xtix, 'XTickLabel',xtix*1.5)
                    xticks([0,33,73,105])
                    xticklabels({'0','50','110','157'})
                    xlabel('Position (cm)');
                    set(gca,'FontSize',22);
                    
                    
                    1;
                end
                
                % Check that there is still a peak at the manipulated landmark position
                valid_peak_positions = [manipulated_landmark_position-14:manipulated_landmark_position+13];
                if sum(ismember(peaks,valid_peak_positions))
                    
                    % Check that the landmark response is persistent.
                    rastermap_at_landmark = rmaps2(:,valid_peak_positions,r);
                    lap_responses_max = nanmax(rastermap_at_landmark'); % Each element in array is the max response at landmark at a given lap
                    lap_responses_mean = nanmean(rastermap_at_landmark'); % Each element in array is the max response at landmark at a given lap
                    
                    % At least 20% of the laps there needs to be a response
                    if sum(lap_responses_max > 0)/length(lap_responses_max) > 0.2
                        % Add cell to trace cells
                        tcs = [tcs,r];
                    end
                end
                
            end
            
            
            % What cells are both LCS and TCS
            tcs = tcs(ismember(tcs,lcs));
            
            
        case 'new'
            
            
            tcs = [];
            rmaps2 = rmaps2(:,:,lcs);
            
            for r = 1:size(rmaps2,3)
                
                % Create lap-averaged map of roi
                avg_map = nanmean(rmaps2(:,:,r));
                
                % Find bins with least activity
                [~,sort_inds] = sort(avg_map);
                
                rmap = rmaps2(:,:,r);
                
                baseline_rmap = rmap(:,sort_inds(1:10));
                baseline_mean = nanmean(baseline_rmap(:));
                baseline_std = nanstd(baseline_rmap(:));
                
                threshold = baseline_mean + (baseline_std*3);
                
                % The lc is considered a tc if it is above threshold in 10%
                % of laps.
                window_of_activity = manipulated_landmark_position-10:manipulated_landmark_position+10;
                
                reliability = 0;
                for lap = 1:size(rmap)
                    if max(rmap(lap,window_of_activity))>threshold
                        reliability = reliability + 1;
                    end
                end
                
                reliability = reliability / size(rmap,1);
                
                if reliability > 0.2
                    tcs = [tcs,r];
                end
            end
            
            % Get tcs back on the correct indexing
            tcs = lcs(tcs);
            
            % Set lcs to be cells that are not tcs
            lcs = lcs(ismember(lcs,tcs) == 0);
    end
    
    
    
else
    
    lcs = [];
    tcs = [];
    
end

end