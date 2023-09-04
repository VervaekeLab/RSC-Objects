function placeCellPeakDistribution(db,varargin)
% placeCellPeakDistribution
% INPUT
%   db: Database containing all sessions that you want to be quantified.

% Params
params.raster_signal = 'deconv';

peaks = [];

for session = 1:length(db)
    
    % Get place cells for session
    place_cells = db(session).place_cells;
    
    % Grab rastermaps for this session
    rmaps = db(session).roimaps.(params.raster_signal);
   
    % Grab landmark on laps
    landmarks_on = db(session).motor_on_laps;
    
    % Use only laps where the motor cues are on and Use only place cell rastermaps
    rmaps = rmaps(landmarks_on,:,place_cells);
    
    % For each roi in rastermap do analysis
    for r = 1:size(rmaps,3)
        
%         figure(221);
%         clf;
%         subplot(2,1,1);
%         imagesc(rmaps(:,:,r));
%         colormap('jet');
%         subplot(2,1,2);
%         plot(nanmean(rmaps(:,:,r)));
%         xlim([1,105])

        % Estimate place cell statistics
        [~,stats] = sb.classify.placeTunedRoisSimple(rmaps(:,:,r),'min_peak_distance',10);
        
        peaks = [peaks,stats.field_start{:}];
    end


%     for x = 1:40
%        subplot(4,10,x);
%        sb.plot.rasterMaps(rmaps(:,:,1:40));
%        colormap('jet');
%     end



end

figure(12912);
clf;
histogram(peaks,105);


end