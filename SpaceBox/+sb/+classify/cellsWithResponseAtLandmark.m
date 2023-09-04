function ccs = cellsWithResponseAtLandmark(rmaps, varargin)
%% Find cells that have a response at either of the landmark positions.

%% Parameters
params.visualize = false;
params.max_peak_width = 50;
params.min_peak_dist = 20;
params.n_peak = 4;
params.min_peak_width = 4;
params.high_low_difference = 0.25;
params.reliability = 0.2;

% Update parameters
params = sb.helper.updateParameterStruct(params,varargin);

%% Find landmark cells based on rmaps
% Init
ccs_stats = struct();
ccs = [];

for r = 1:size(rmaps,3)
    
    % Create lap-averaged map of roi
    avg_map = nanmean(rmaps(:,:,r));
    
    % Smooth the averaged map
    avg_map_smoothed = smoothdata(avg_map,'Gaussian',20);
    
    % Set parameters for find peaks
    % Set a minimum threshold for a peak to be included
    peak_threshold = min(avg_map_smoothed)+(params.high_low_difference*(max(avg_map_smoothed)-min(avg_map_smoothed))); % Threshold is set as 20 % of the difference between the highest and lowest activity of the position activity map.

    %[~,sort_ind] = sort(avg_map_smoothed);
    %peak_threshold = nanmean(avg_map_smoothed(sort_ind(1:10))) + (nanstd(avg_map_smoothed(sort_ind(1:10)))*3);
    [~,peaks] = findpeaks(avg_map_smoothed,'NPeaks',params.n_peak,'MinPeakHeight',peak_threshold,'MinPeakDistance',params.min_peak_dist,'SortStr','descend','MaxPeakWidth',params.max_peak_width,'MinPeakWidth',params.min_peak_width);
    
    % Find the two peaks close to the landmarks
    if length(peaks)>2
        peaks = [peaks(ismember(peaks,23:43)),peaks(ismember(peaks,63:83))];
    end
    
    if params.visualize
        figure(1); clf;
        subplot(1,4,1:3);
        findpeaks(avg_map_smoothed,'NPeaks',params.n_peak,'MinPeakHeight',peak_threshold,'MinPeakDistance',params.min_peak_dist,'SortStr','descend','MaxPeakWidth',params.max_peak_width,'MinPeakWidth',params.min_peak_width);
        hold on;
        
        line([33,33],[-1000,1000]);
        line([73,73],[-1000,1000]);
        line([93,93],[-1000,1000]);
        line([0,105],[peak_threshold,peak_threshold],'Color','r','LineStyle','-','LineWidth',2);
        plot(avg_map,'Color','k','LineStyle','--','LineWidth',2);
        
        xtix = get(gca, 'XTick');
        set(gca, 'XTick',xtix, 'XTickLabel',xtix*1.5)
        xticks([0,33,73,105])
        xticklabels({'0','50','110','157'})
        xlabel('Position (cm)');
        set(gca,'FontSize',22);
        
        subplot(1,4,4);        
        imagesc(rmaps(:,:,r));
        title('NO');
        1;
        xtix = get(gca, 'XTick');
        set(gca, 'XTick',xtix, 'XTickLabel',xtix*1.5)
        xticks([0,33,73,105])
        xticklabels({'0','50','110','157'})
        xlabel('Position (cm)');
        set(gca,'FontSize',22);
        1;
    end
    
    % Only include cells with two peaks
    if length(peaks) > 1
        
        if (length(peaks(ismember(peaks,23:43)))) || (length(peaks(ismember(peaks,63:83))))
        
            % Compute reliability at each positions
            ccs_stats(r).reliability_first  = sum(max(rmaps(:,23:43,r)')>0) / size(rmaps,1);        
            ccs_stats(r).reliability_second = sum(max(rmaps(:,63:83,r)')>0) / size(rmaps,1);
            
            % Check that the reliability is big enough on either of the two
            if ccs_stats(r).reliability_first > params.reliability
                css_stats(r).response_at_first = 1;
                
            else
                css_stats(r).response_at_first = 0;
            end
            
            if ccs_stats(r).reliability_second > params.reliability
                css_stats(r).response_at_second = 1;
            else
                css_stats(r).response_at_second = 0;
            end
            
            % If the cell responds reliably at either position add it to
            % list of conjunctive cells.
            if css_stats(r).response_at_first || css_stats(r).response_at_second
               ccs = [ccs, r]; 
            end
            
        end
          
    end
    
end

end