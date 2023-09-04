% Show a pie chart of the percentage of persistent vs non-persistent
% landmark cells

figure(3);
p = pie([sum(data.n_tcs),sum(data.n_lcs)]);
p(1).FaceColor = 'k';
p(3).FaceColor = [0.9,0.9,0.9];


%% Show an overlap between predictive cells and persitency
predictive_cells = predictiveCells(mData);

n_predictive_cells = length(predictive_cells);
n_persistent_cells = sum(data.n_tcs);
n_non_persistent_cells = sum(data.n_lcs);

% Check how many overlap between categories
n_persistent_and_non_persistent_cells = 0; % by definition it will be zero 
n_predictive_and_persistent = 0;
n_predictive_and_non_persistent = 0;
n_persistent_and_non_persistent_cells_and_predictive = 0;

for c = 1:length(predictive_cells.session)
   n_predictive_and_persistent = n_predictive_and_persistent + ismember(predictive_cells.ind(c),mData(predictive_cells.session(c)).rmaps.tcs.ind);
   n_predictive_and_non_persistent = n_predictive_and_non_persistent + ismember(predictive_cells.ind(c),mData(predictive_cells.session(c)).rmaps.lcs.ind);
end

figure(4); clf;

subplot(2,1,1);
pie([n_predictive_and_persistent,n_persistent_cells])
subplot(2,1,2);
pie([n_predictive_and_non_persistent,n_non_persistent_cells])


%% Print out the percentage of predictive cells that are peristent and non-persistent
pred_persistent_ratio = (n_predictive_and_persistent/(n_predictive_and_persistent+n_predictive_and_non_persistent))*100;
pred_non_persistent_ratio = (n_predictive_and_non_persistent/(n_predictive_and_persistent+n_predictive_and_non_persistent))*100;
fprintf('Predictive cells that are persistent: %.1f %%\n',pred_persistent_ratio)
fprintf('Predictive cells that are non-persistent: %.1f %%\n',pred_non_persistent_ratio)


%% Print the percentage of non-predictive cells that are persistent and non-persistent

% For each session, grab all cells that are non predictive
n_nonpredictive_and_non_persistent = 0;
n_nonpredictive_and_persistent = 0;

for s = 1:length(mData)
    
    all_cells = [mData(s).rmaps.lcs.ind,mData(s).rmaps.tcs.ind];
    pred_cells = predictive_cells.ind(predictive_cells.session == s);
    non_pred_cells = all_cells(ismember(all_cells, pred_cells) == 0);
    
    n_nonpredictive_and_non_persistent = n_nonpredictive_and_non_persistent + sum(ismember(mData(s).rmaps.lcs.ind,non_pred_cells));
    n_nonpredictive_and_persistent = n_nonpredictive_and_persistent + sum(ismember(mData(s).rmaps.tcs.ind,non_pred_cells));

end

non_pred_persistent_ratio = (n_nonpredictive_and_persistent/(n_nonpredictive_and_persistent+n_nonpredictive_and_non_persistent))*100;
non_pred_non_persistent_ratio = (n_nonpredictive_and_non_persistent/(n_nonpredictive_and_non_persistent+n_nonpredictive_and_persistent))*100;

fprintf('Non-predictive cells that are persistent: %.1f %%\n',non_pred_persistent_ratio)
fprintf('Non-predictive cells that are non-persistent: %.1f %%\n',non_pred_non_persistent_ratio)


%% Plot distribution of predictive and persistency
figure(10);clf;
b = barh([non_pred_persistent_ratio,non_pred_non_persistent_ratio;pred_persistent_ratio,pred_non_persistent_ratio],'stacked');
b(1).FaceColor = [0.2,0.2,0.2];
b(2).FaceColor = [0.8,0.8,0.8];

% Set figure size
set(gcf,'renderer', 'painters', 'Position', [-1200,1000,650,300])


%% Helper functions
function predictive_cells = predictiveCells(mData)
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
    predictive_cells_ind = indices(first_peak_cells);
    predictive_cells_session = [predictive_cells_session, session(second_peak_cells)];
    predictive_cells_ind = [predictive_cells_ind, indices(second_peak_cells)];
    
    predictive_cells.session = predictive_cells_session;
    predictive_cells.ind = predictive_cells_ind;
    

end