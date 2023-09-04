

% find all predictive cells in jitter
amap_object_first = [];
amap_object_second = [];
amap_object_third = [];
amap_object_fourth = [];
no_all_cells = 0;
for s = 1:length(mData_jitter)
     amap_object_first = cat(1,amap_object_first,sb.generate.averagedTuningFromRasterMaps(mData_jitter(s).sData.rmaps.deconv_firstpair(:,:,mData_jitter(s).rmaps_first)));
     amap_object_second = cat(1,amap_object_second,sb.generate.averagedTuningFromRasterMaps(mData_jitter(s).sData.rmaps.deconv_secondpair(:,:,mData_jitter(s).rmaps_second)));
     amap_object_third = cat(1,amap_object_third,sb.generate.averagedTuningFromRasterMaps(mData_jitter(s).sData.rmaps.deconv_thirdpair(:,:,mData_jitter(s).rmaps_third)));
     amap_object_fourth = cat(1,amap_object_fourth,sb.generate.averagedTuningFromRasterMaps(mData_jitter(s).sData.rmaps.deconv_fourthpair(:,:,mData_jitter(s).rmaps_fourth)));

    no_all_cells = no_all_cells+ size(mData_jitter(s).sData.rmaps.deconv,3);
end


amap_object_first  = sb.sort.rasterMapByPeak(amap_object_first ,'normalize_type','zscore');
amap_object_second  = sb.sort.rasterMapByPeak(amap_object_second ,'normalize_type','zscore');
amap_object_third  = sb.sort.rasterMapByPeak(amap_object_third ,'normalize_type','zscore');
amap_object_fourth  = sb.sort.rasterMapByPeak(amap_object_fourth ,'normalize_type','zscore');

first_landmark_position = 22; second_landmark_position_threshold = round(55-0.5*(55-22)); 
predictive_rmaps_first  = plot_predictive(gca,amap_object_first,first_landmark_position,second_landmark_position_threshold);

first_landmark_position = 32; % This is set to be the position (in bins) where the motor first starts to move. Whisker touches don't happen until at least bin 32.
second_landmark_position_threshold = round(65-0.5*(65-32));
predictive_rmaps_second  = plot_predictive(gca,amap_object_second,first_landmark_position,second_landmark_position_threshold);

first_landmark_position = 41; % This is set to be the position (in bins) where the motor first starts to move. Whisker touches don't happen until at least bin 32.
second_landmark_position_threshold = round(75-0.5*(75-41)); % This is set to be the position (in bins) where the motor first starts to move. Whisker touches don't happen until at least bin 72.
predictive_rmaps_third  = plot_predictive(gca,amap_object_third,first_landmark_position,second_landmark_position_threshold);

first_landmark_position = 52; % This is set to be the position (in bins) where the motor first starts to move. Whisker touches don't happen until at least bin 32.
second_landmark_position_threshold = round(52-0.5*(85-52)); % This is set to be the position (in bins) where the motor first starts to move. Whisker touches don't happen until at least bin 72.
predictive_rmaps_fourth  = plot_predictive(gca,amap_object_fourth,first_landmark_position,second_landmark_position_threshold);

predictive_jitter_prct = 100*((size(predictive_rmaps_first,1)+size(predictive_rmaps_second,1)...
    +size(predictive_rmaps_third,1)+size(predictive_rmaps_fourth,1))/no_all_cells);

% familiar
predicitve_familiar = n_predictive_cells ./ sum([n_predictive_cells,n_nonpredictive_cells]) * 100;

figure()
bar([predicitve_familiar,predictive_rmaps_introduction_prct,predictive_jitter_prct],'FaceColor',[0 0 0])
xticks([1:3])
xticklabels({'familiar', 'introduction','jitter'})
set(gca,'FontSize',20,'FontName','Arial');
ylabel('Predictive cells (%)')
ylim([0 30])
yticks([0:10:30])
box off

function predictive_rmaps  = plot_predictive(gca,amap,first_landmark_position,second_landmark_position)

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

end
