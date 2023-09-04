


% [~,motors] = split_rmaps_jitter(mData_jitter(1).sData);

pos = mData_jitter(1).sData.behavior.wheelPosDsBinned;
lap = mData_jitter(1).sData.behavior.wheelLapDsBinned;
run = mData_jitter(1).sData.behavior.runSpeedDs;

right_motor = mData_jitter(1).sData.daqdata.right_motor_speed(mData_jitter(1).sData.daqdata.frameIndex);
rmaps.motor_right = sb.generate.rasterMaps(pos,lap,right_motor,'smoothing',0);

trialType = ones(size(rmaps.motor_right));
trialType(rmaps.motor_right>0) = 0;
trialTypeinvert= ones(size(trialType));
for i = 2:size(trialType,1)-2
    trialType(i,:)= [0,diff(trialType(i,:))];
    idx = find(trialType(i,:) > 0);
    idx_all= [idx(1)-1:idx(1)+2, idx(2)-1:idx(2)+2];
    trialTypeinvert(i,idx_all) = 0;
end

figure()
imagesc(trialTypeinvert(2:13,:))
ylabel('Lap')
set(gca,'FontName','Arial', 'FontSize',14)
xticks([1,33,73,105])
xticklabels({'0','50','110','157'})
colormap('gray'); 
yticks([1 12])


figure()

scatter(ones(1,6)',100*data_jitter.prct_lcs_in_all',50,'k')
hold on
scatter(1,100*nanmean(data_jitter.prct_lcs_in_all'),70,'k','filled')
errorbar(1,100*nanmean(data_jitter.prct_lcs_in_all'),nanstd(100*data_jitter.prct_lcs_in_all')/sqrt(length(data_jitter.prct_lcs_in_all)),'k','LineWidth',1.5)
ylabel('Object resposinve cells in all 4 (%)')
set(gca,'FontName','Arial', 'FontSize',14)
xticks([])
yticks([0:10:50])


amap_object_responsive_in_all_pos1 = [];
amap_object_responsive_in_all_pos2 = [];
amap_object_responsive_in_all_pos3 = [];
amap_object_responsive_in_all_pos4 = [];
for s = 1:length(mData_jitter)

     amap_object_responsive_in_all_pos1 = cat(1,amap_object_responsive_in_all_pos1,sb.generate.averagedTuningFromRasterMaps(mData_jitter(s).sData.rmaps.deconv_firstpair(:,:,data_jitter.lcs_in_all{s})));
     amap_object_responsive_in_all_pos2 = cat(1,amap_object_responsive_in_all_pos2,sb.generate.averagedTuningFromRasterMaps(mData_jitter(s).sData.rmaps.deconv_secondpair(:,:,data_jitter.lcs_in_all{s})));
     amap_object_responsive_in_all_pos3 = cat(1,amap_object_responsive_in_all_pos3,sb.generate.averagedTuningFromRasterMaps(mData_jitter(s).sData.rmaps.deconv_thirdpair(:,:,data_jitter.lcs_in_all{s})));
     amap_object_responsive_in_all_pos4 = cat(1,amap_object_responsive_in_all_pos4,sb.generate.averagedTuningFromRasterMaps(mData_jitter(s).sData.rmaps.deconv_fourthpair(:,:,data_jitter.lcs_in_all{s})));
    
end

amap_object_responsive_in_all_pos1  = sb.sort.rasterMapByPeak(amap_object_responsive_in_all_pos1 ,'normalize_type','zscore');
amap_object_responsive_in_all_pos2  = sb.sort.rasterMapByPeak(amap_object_responsive_in_all_pos2 ,'normalize_type','zscore');
amap_object_responsive_in_all_pos3  = sb.sort.rasterMapByPeak(amap_object_responsive_in_all_pos3 ,'normalize_type','zscore');
amap_object_responsive_in_all_pos4  = sb.sort.rasterMapByPeak(amap_object_responsive_in_all_pos4 ,'normalize_type','zscore');

figure()
subplot(4,1,1)
imagesc(amap_object_responsive_in_all_pos1,[0,max_z_value])
xticks([0,33,73,105])
xticklabels({'0','50','110','157'})
set(gca,'FontSize',20);
line([22,22],[0,2000],'Color','w');
line([55,55],[0,2000],'Color','w');
ylabel('Cell #')


first_landmark_position = 22; % This is set to be the position (in bins) where the motor first starts to move. Whisker touches don't happen until at least bin 32.
second_landmark_position_threshold = round(55-0.5*(55-22)); % This is set to be the position (in bins) where the motor first starts to move. Whisker touches don't happen until at least bin 72.

[predictive_rmaps_pos1, predictive_prct_pos1]= plot_predictive(gca,amap_object_responsive_in_all_pos1,first_landmark_position ,second_landmark_position_threshold,data_jitter);
% Landmark onset
% Set correct size of the figure


subplot(4,1,2)
imagesc(amap_object_responsive_in_all_pos2,[0,max_z_value])
xticks([1,33,73,105])
xticklabels({'0','50','110','157'})
set(gca,'FontSize',20);
line([32,32],[0,2000],'Color','w');
line([65,65],[0,2000],'Color','w');

first_landmark_position = 32; % This is set to be the position (in bins) where the motor first starts to move. Whisker touches don't happen until at least bin 32.
second_landmark_position_threshold = round(65-0.5*(65-32)); % This is set to be the position (in bins) where the motor first starts to move. Whisker touches don't happen until at least bin 72.

[predictive_rmaps_pos2,predictive_prct_pos2] = plot_predictive(gca,amap_object_responsive_in_all_pos2,first_landmark_position ,second_landmark_position_threshold,data_jitter);


subplot(4,1,3)
imagesc(amap_object_responsive_in_all_pos3,[0,max_z_value])
xticks([1,33,73,105])
xticklabels({'0','50','110','157'})
set(gca,'FontSize',20);
line([42,42],[0,2000],'Color','w');
line([75,75],[0,2000],'Color','w');

first_landmark_position = 41; % This is set to be the position (in bins) where the motor first starts to move. Whisker touches don't happen until at least bin 32.
second_landmark_position_threshold = round(75-0.5*(75-41)); % This is set to be the position (in bins) where the motor first starts to move. Whisker touches don't happen until at least bin 72.

[predictive_rmaps_pos3 ,predictive_prct_pos3]= plot_predictive(gca,amap_object_responsive_in_all_pos3,first_landmark_position ,second_landmark_position_threshold,data_jitter);


subplot(4,1,4)
imagesc(amap_object_responsive_in_all_pos4,[0,max_z_value])
xticks([1,33,73,105])
xticklabels({'0','50','110','157'})
set(gca,'FontSize',20);
line([52,52],[0,2000],'Color','w');
line([85,85],[0,2000],'Color','w');


first_landmark_position = 52; % This is set to be the position (in bins) where the motor first starts to move. Whisker touches don't happen until at least bin 32.
second_landmark_position_threshold = round(52-0.5*(85-52)); % This is set to be the position (in bins) where the motor first starts to move. Whisker touches don't happen until at least bin 72.

[predictive_rmaps_pos4 ,predictive_prct_pos4]= plot_predictive(gca,amap_object_responsive_in_all_pos4,first_landmark_position ,second_landmark_position_threshold,data_jitter);

precitive_prct_jitter = nansum([predictive_prct_pos1;predictive_prct_pos2;predictive_prct_pos3;predictive_prct_pos4],1);

nanmean(precitive_prct_jitter)
nanstd(precitive_prct_jitter)/sqrt(6)


function [predictive_rmaps,predictive_prct]  = plot_predictive(gca,amap,first_landmark_position,second_landmark_position,data_jitter)

    [~,first_peaks] = max(amap');

    % Sort session and indices based on the sorted map

    first_peak_cells = first_peaks<first_landmark_position & first_peaks>16;
    first_peak_cells = find(first_peak_cells);

    % Find first cell where we are within the second landmark regime
    first_cell_second_landmark = find(first_peaks>second_landmark_position);
    first_cell_second_landmark = first_cell_second_landmark(1);

    [~,second_peaks] = max(amap(first_cell_second_landmark:end,:)');

    second_peak_cells = second_peaks<second_landmark_position;
    second_peak_cells = find(second_peak_cells)+first_cell_second_landmark;
    %sb.plot.rasterMaps(amap_lcs_sorted([first_peak_cells,second_peak_cells],:))

    % Save rmaps for 50 laps of all predictive cells 
    % predictive_rmaps = zeros(50,105);
    predictive_rmaps= [];
    for c = 1:length(first_peak_cells)
       predictive_rmaps = [predictive_rmaps; amap(first_peak_cells(c),:)];
    end

    for c = 1:length(second_peak_cells)
       predictive_rmaps = [predictive_rmaps; amap(second_peak_cells(c),:)];
    end

    line_width = 1;
    line_color = 'r';

    % Top
    if ~isempty(first_peak_cells)
        line(gca,[0,105],[first_peak_cells(1),first_peak_cells(1)],'Color',line_color,'LineWidth',line_width);
        line(gca,[0,105],[first_peak_cells(end),first_peak_cells(end)],'Color',line_color,'LineWidth',line_width);
        line(gca,[1,1],[first_peak_cells(1),first_peak_cells(end)],'Color',line_color,'LineWidth',line_width);
        line(gca,[105,105],[first_peak_cells(1),first_peak_cells(end)],'Color',line_color,'LineWidth',line_width);
    end
    % Bottom
    if ~isempty(second_peak_cells)
        line(gca,[0,105],[second_peak_cells(1),second_peak_cells(1)],'Color',line_color,'LineWidth',line_width);
        line(gca,[0,105],[second_peak_cells(end),second_peak_cells(end)],'Color',line_color,'LineWidth',line_width);
        line(gca,[1,1],[second_peak_cells(1),second_peak_cells(end)],'Color',line_color,'LineWidth',line_width);
        line(gca,[105,105],[second_peak_cells(1),second_peak_cells(end)],'Color',line_color,'LineWidth',line_width);
    end
    
    session = [];
    for s = 1:length(data_jitter)
        session = [session, ones(1,length(data_jitter.lcs_in_all{s}))*s];
    end
    
    predictive_cells_session = session(first_peak_cells);
    predictive_cells_session = [predictive_cells_session, session(second_peak_cells)];
    
    if ~isempty(predictive_cells_session)
        for s = 1:length(data_jitter)
            predictive_prct(s) = 100*find(predictive_cells_session ==s)/length(data_jitter.n_lcs_in_all);
        end
    else
        predictive_prct = zeros(1,length(data_jitter.n_lcs_in_all));
    end
  

end