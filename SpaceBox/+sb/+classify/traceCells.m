function tcs = traceCells(motor_on_rastermaps,motor_off_rastermaps)
%% traceCells
% Find all trace cells comparing the two 3D rastersmaps given as input.
%
% INPUT
%   -- Required
%       motor_on_rastermaps:
%       motor_off_rastermaps:
%
% OUTPUT
%

params.visualize = false;

%% Find trace cells
% Find place cells
pcs = find(sb.classify.placeTunedRoisSimple(motor_on_rastermaps,'reliability_criteria',0.33));

tcs = [];
for i = 1:length(pcs)
    
    % Create a smoothed averaged activity profile
    motor_on_map = nanmean(smoothdata(motor_on_rastermaps(:,:,pcs(i)),2,'gaussian',20));
    motor_off_map = nanmean(smoothdata(motor_off_rastermaps(:,:,pcs(i)),2,'gaussian',20));
    
    % Calculate minimum threshold for peaks to be detected based on Mao criteria
    pf_threshold = min(motor_on_map)+(0.3*(max(motor_on_map)-min(motor_on_map))); % Threshold is set as n % of the difference between the highest and lowest activity of the position activity map.

    % Find peaks
    [~,peaks_on] = findpeaks(motor_on_map,'NPeaks',3, ...
        'MinPeakDistance', 12, ...
        'SortStr','descend', ...
        'MinPeakHeight',pf_threshold);
    
    
    % Check that peaks are not at edge of track
    % ... 
    
    % Visualize
    if params.visualize
        
        figure(1);clf;
        plot(motor_on_map,'Color','k','LineWidth',3);
        hold on
        plot(motor_off_map,'Color','r','LineWidth',3);
        hold off
        
        figure(2); clf;
        findpeaks(motor_on_map,'NPeaks',3, ...
            'MinPeakDistance', 12, ...
            'SortStr','descend', ...
            'MinPeakHeight',pf_threshold);
        
    end
    
    % Check that there are more than one peak in motor on condition
    if length(peaks_on) > 1
        peaks_on = peaks_on(1:2);
        
        % Is the motor off peaks at the same position?
        [~,peaks_off] = findpeaks(motor_off_map,'NPeaks',3, ...
            'MinPeakDistance', 12, ...
            'SortStr','descend');
        
        % Are there more than 1 peak in the motor off condition
        if length(peaks_off) > 1
            peaks_off = peaks_off(1:2);
            
            % Set valid peak positions to be +/- 5 bins from the position
            % in the motor on condition
            valid_peaks = 0;
            valid_pos = [];
            for p = 1:2
                valid_pos = [valid_pos,peaks_on(p)-5:peaks_on(p)+5];
            end
            
            if sum(peaks_off(1)==valid_pos)
                valid_peaks = valid_peaks + 1;
            end
            if sum(peaks_off(2)==valid_pos)
                valid_peaks = valid_peaks + 1;
            end
            
            if valid_peaks == 2
                tcs = [tcs,i];
                
            else
                % if there are not two valid peaks check if the peaks are 40 ish
                % bins away from each other
                try
                    on_diff = max(peaks_on)-min(peaks_on);
                    off_diff = max(peaks_off)-min(peaks_off);
                    
                    if (on_diff > 33) & (on_diff < 47)
                        
                       if  (off_diff > 33) & (off_diff < 47)
                          tcs = [tcs,i]; 
                       end
                    end 
                end
            end
        end
    end
end





end