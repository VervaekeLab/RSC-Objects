function lcs = landmarkCells(rmaps)

params.visualize = false;

% Landmarks are at 33 and 73 bin. I take 20 cm before and after which is
% 20/1.5 = 13.3 bins
valid_landmark_positions = [33-13:33+13,73-13:73+13];

average_rmaps = sb.generate.averagedTuningFromRasterMaps(rmaps,'smoothing',20);
av_rmaps = sb.generate.averagedTuningFromRasterMaps(rmaps,'smoothing',15);
lcs = [];

n_landmark_peaks = [];
for c = 1:size(rmaps,3)
    
    if c == 218
        %params.visualize = true;
    end
    
    if sb.classify.placeTunedRoisSimple(rmaps(:,:,c),'reliability_criteria',0.3)
        peak_threshold = min(average_rmaps(c,:))+(0.25*(max(average_rmaps(c,:))-min(average_rmaps(c,:)))); % Threshold is set as n % of the difference between the highest and lowest activity of the position activity map.
        
        [~,peaks] = findpeaks(average_rmaps(c,:),'NPeaks',3,'MinPeakHeight',peak_threshold,'MinPeakDistance',25);
        
        
        n_landmark_peaks(c) = 0;
        
        if sum(ismember(peaks,valid_landmark_positions)) == 2
            n_landmark_peaks(c) = 2;
        end
        
        
        peaks = peaks(peaks<88);
        peaks = peaks(peaks>13);
        
        if length(peaks) == 2
            
            % this is a landmark cell
            lcs = [lcs,c];
            
            if params.visualize
                figure(1);
                
                clf;
                findpeaks(average_rmaps(c,:),'NPeaks',3,'MinPeakHeight',peak_threshold,'MinPeakDistance',25);
                hold on;
                line([33,33],[-100,1000]);
                line([73,73],[-100,1000]);
                line([93,93],[-100,1000]);
                hold off;
                
                
                figure(2); clf;
                imagesc(rmaps(:,:,c));
                
                hold on;
                line([33,33],[-100,1000]);
                line([73,73],[-100,1000]);
                line([93,93],[-100,1000]);
                hold off;
                fprintf('%i: %i cm\n',c,round((max(peaks)-min(peaks))*1.5));
                
                
            end
            
        end
    end
end






end