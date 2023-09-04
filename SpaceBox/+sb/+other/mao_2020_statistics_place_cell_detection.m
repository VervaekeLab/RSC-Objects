function pcs = mao_2020_statistics_place_cell_detection(sData)
% Method used by Mao (2020, Cell) for detecting place tuned cells. In their
% methods section they write that they tried two different methods. One was
% the following shuffling method which is implemented, and the other was
% using spatial information. They state that the methods gave similar
% results and so they choose the shuffling method.
%
% Note Andreas 03.2020: I tried the method on my dataset and it detects
% place cells but also cells that I would not considered place tuned by
% looking at the responses.

% Get data
pos = sData.behavior.wheelPosDsBinned;
lap = sData.behavior.wheelLapDsBinned;
deconv = sData.imdata.roiSignals(2).deconv; % Deconvolved signal using NNMF (Pnevmatikakis)


% --- This part is skipped, only necessary for my recordings where I do
% % manipulations halfway through.
% 
% % Cut recording when the first landmark is omitted for the first time
% % ...
% lap_bigger_than_50 = find(lap>50);
% lap_bigger_than_50 = lap_bigger_than_50(1);
% 
% pos = pos(1:lap_bigger_than_50);
% lap = lap(1:lap_bigger_than_50);
% 
% new_deconv = [];
% for r = 1:size(deconv,1)
%     new_deconv(r,:) = deconv(r,1:lap_bigger_than_50);
% end
% 
% deconv = new_deconv;
% ----

% Create occupancy-normalized activity maps
rmaps = sb.generate.rasterMaps(pos,lap,deconv,'smoothing',3); % Smoothing 3 bins

% For each roi, do shuffling and create new rastermaps and compare with
% original rmap
pcs = [];
for roi = 1:size(rmaps,3)
    
    % Shuffle 1000 times
    shuffled_distribution = [];
    for shuffle = 1:1000
        % We create a new rastermap for the roi by circshifting the deconv
        % signal
        shift_size = randi((size(deconv,2)-(40*31)),1)+ (20*31); % shift size is between 20 and session duration minus 20 s.
        roi_deconv = circshift(deconv(roi,:),shift_size);
        shifted_rmap = sb.generate.rasterMaps(pos,lap,roi_deconv,'smoothing',3);
        shuffled_distribution(shuffle) = nanmean(nanmean(shifted_rmap));
    end
    
    % Calculate 97.5 percentile of the shuffled distribution
    upper_pct = prctile(shuffled_distribution,97.5);
    
    % According to Mao, a cell is classified as spatially tuned if its
    % lower bound (mean - SEM) in any position is bigger than the the 97.5
    % percentile of the shuffled distribution.
    
    avg_map = nanmean(rmaps(:,:,roi));
    avg_map_mean = mean(avg_map);
    avg_map_sem = std(avg_map)/sqrt(length(avg_map));
    
    % Add roi to list of place cells if criteria is fullfilled.
    if (avg_map_mean - avg_map_sem) > upper_pct
       pcs = [pcs, roi];
    end
    
end


end