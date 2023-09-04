% Load the main dataset
clear;
[mData,data] = sb.database.load.rsc_manipulation;

%% Get average tuning of all landmark cells and find predictive cells
% Parameters
first_landmark_position = 32; % This is set to be the position (in bins) where the motor first starts to move. Whisker touches don't happen until at least bin 32.
second_landmark_position = 72; % This is set to be the position (in bins) where the motor first starts to move. Whisker touches don't happen until at least bin 72.

amap_lcs = [];
signal_type = "deconv_motor_on";
signal_type_visualization = "deconv_motor_on";
session = [];
indices = [];
for s = 1:length(mData)
    inds = [mData(s).rmaps.lcs.ind,mData(s).rmaps.tcs.ind];
    session = [session, ones(1,length(inds))*s];
    indices = [indices,inds];
    if s == 1
        amap_lcs = sb.generate.averagedTuningFromRasterMaps(cat(3,mData(s).rmaps.lcs.(signal_type),mData(s).rmaps.tcs.(signal_type)));    
    else
        amap_lcs = cat(1,amap_lcs,sb.generate.averagedTuningFromRasterMaps(cat(3,mData(s).rmaps.lcs.(signal_type),mData(s).rmaps.tcs.(signal_type))));
    end
    
end

[amap_lcs_sorted,sort_ind] = sb.sort.rasterMapByPeak(amap_lcs,'normalize_type','zscore');
[~,first_peaks] = max(amap_lcs_sorted');

% Sort session and indices based on the sorted map
session = session(sort_ind);
indices = indices(sort_ind);

first_peak_cells = first_peaks<first_landmark_position;
first_peak_cells = find(first_peak_cells);

% Find first cell where we are within the second landmark regime
first_cell_second_landmark = find(first_peaks>50);
first_cell_second_landmark = first_cell_second_landmark(1);

[~,second_peaks] = max(amap_lcs_sorted(first_cell_second_landmark:end,:)');

second_peak_cells = second_peaks<second_landmark_position;
second_peak_cells = find(second_peak_cells)+first_cell_second_landmark;
%sb.plot.rasterMaps(amap_lcs_sorted([first_peak_cells,second_peak_cells],:))

% Save rmaps for 50 laps of all predictive cells 
predictive_rmaps = zeros(50,105);

for c = 1:length(first_peak_cells)
   predictive_rmaps = cat(3,predictive_rmaps, mData(session(c)).rmaps.(signal_type_visualization)(1:50,:,indices(c)));
end

for c = 1:length(second_peak_cells)
   predictive_rmaps = cat(3,predictive_rmaps, mData(session(second_peak_cells(c))).rmaps.(signal_type_visualization)(1:50,:,indices(second_peak_cells(c))));  
end

predictive_rmaps = predictive_rmaps(:,:,2:end);

predictive_cells_session = session(first_peak_cells);
predictive_cells_inds = indices(first_peak_cells);
predictive_cells_session = [predictive_cells_session, session(second_peak_cells)];
predictive_cells_inds = [predictive_cells_inds, indices(second_peak_cells)];

for i = 1:length(mData)
    predictive_cells_familiar(i) = 100*length(find(predictive_cells_session == i))/length(unique([mData(i).rmaps.lcs.ind,mData(i).rmaps.tcs.ind]));
end

nanmean(predictive_cells_familiar)
nanstd(predictive_cells_familiar)/sqrt(length(mData))
%% Plot all predictive cells
%sb.plot.rasterMaps(predictive_rmaps,'cmap_max_type','peak','add_lines_vertical',[30,70])