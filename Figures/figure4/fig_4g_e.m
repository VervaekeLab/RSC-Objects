


landmark_width = 3;

landmarks = mData_introduction(2).sData.landmarks.motor.on(:,1);
landmarks = landmarks(find(~isnan(landmarks)));
sorted_landmarks = flipud(sort(landmarks));

% Set a value for each lap where a landmark is present
motormap = zeros(105,length(landmarks));
% motormap(73:73+landmark_width-1,:) = ones(landmark_width,length(landmarks));
motormap_sorted = motormap;
% motormap_sorted(33:33+landmark_width-1,1:sum(landmarks)) = ones(landmark_width,sum(landmarks));

for l = 1:length(landmarks)
   
    if landmarks(l) == 1
        motormap(33:33+landmark_width-1,l) = ones(1,landmark_width);
        motormap(73:73+landmark_width-1,l) = ones(1,landmark_width);

    end
    
end

% Plot
figure();
imagesc(motormap');

cmap = colormap('gray');
colormap(flipud(cmap));
xticks([]);
yticks([1,find(landmarks,1),length(landmarks)])

xticks([1,33,73,105])
xticklabels({'0','50','110','157'})
set(gca,'FontSize',20);
ylabel('Lap #')
xlabel('Position (cm)')


%% Plot all place cells and landmark cells sorted
amap_pcs_on = [];
amap_lcs_on = [];

amap_pcs_off = [];
amap_lcs_off = [];


amap_pcs_on_deconvOff = [];
amap_lcs_on_deconvOff = [];
amap_pcs_off_deconvOn = [];
amap_lcs_off_deconvOn  = [];
max_z_value = 3;

for s = 1:length(mData_introduction)
       session = [session, ones(1,length(mData_introduction(s).rmaps_on.lcs.ind))*s];
%     if s == 1
       amap_pcs_on = cat(1, amap_pcs_on,sb.generate.averagedTuningFromRasterMaps(mData_introduction(s).rmaps_on.pcs.deconv));
       amap_lcs_on = cat(1, amap_lcs_on,sb.generate.averagedTuningFromRasterMaps(mData_introduction(s).rmaps_on.lcs.deconv));
       amap_pcs_off = cat(1, amap_pcs_off,sb.generate.averagedTuningFromRasterMaps(mData_introduction(s).rmaps_off.pcs.deconv));
       
       amap_pcs_on_deconvOff = cat(1,amap_pcs_on_deconvOff ,sb.generate.averagedTuningFromRasterMaps(mData_introduction(s).rmaps_on.pcs.deconv_off));
       amap_lcs_on_deconvOff = cat(1,amap_lcs_on_deconvOff,sb.generate.averagedTuningFromRasterMaps(mData_introduction(s).rmaps_on.lcs.deconv_off));
       

       if ~isempty(mData_introduction(s).rmaps_off.lcs.ind)
            amap_lcs_off = sb.generate.averagedTuningFromRasterMaps(mData_introduction(s).rmaps_off.lcs.deconv);
            amap_lcs_off_deconvOn = cat(1,amap_lcs_off_deconvOn,sb.generate.averagedTuningFromRasterMaps(mData_introduction(s).rmaps_off.lcs.deconv_on));
       end
    
end

session = [];
for s = 1:length(mData_introduction)
       session = [session, ones(1,length(mData_introduction(s).rmaps_on.lcs.ind))*s];
end

% Sort
[amap_lcs_on,sort_idx_on]  = sb.sort.rasterMapByPeak(amap_lcs_on ,'normalize_type','zscore');
amap_lcs_on  = sb.sort.rasterMapByPeak(amap_lcs_on ,'normalize_type','zscore');

amap_lcs_on_deconvOff =normalize(amap_lcs_on_deconvOff(sort_idx_on,:),2,'zscore');

[amap_lcs_off,sort_idx_off]  = sb.sort.rasterMapByPeak(amap_lcs_off ,'normalize_type','zscore');
amap_lcs_off  = sb.sort.rasterMapByPeak(amap_lcs_off ,'normalize_type','zscore');

amap_lcs_off_deconvOn =normalize(amap_lcs_off_deconvOn(sort_idx_off,:),2,'zscore');

% Plot all place and landmark cells sorted

figure(); clf;
subplot(1,2,1);
imagesc(amap_lcs_on_deconvOff,[0,max_z_value]);
hold on
line([32,32],[0,2000],'Color','w');
line([72,72],[0,2000],'Color','w');
yticks([1,size(amap_lcs_on_deconvOff,1)])

xtix = get(gca, 'XTick');
set(gca, 'XTick',xtix, 'XTickLabel',xtix*1.5)
xticks([0,33,73,105])
xticklabels({'0','50','110','157'})
set(gca,'FontSize',20);


subplot(1,2,2)
imagesc(amap_lcs_on,[0,max_z_value]);

hold on
line([32,32],[0,2000],'Color','w');
line([72,72],[0,2000],'Color','w');

xtix = get(gca, 'XTick');
set(gca, 'XTick',xtix, 'XTickLabel',xtix*1.5)
xticks([1,33,73,105])
xticklabels({'0','50','110','157'})
yticks([1,size(amap_lcs_on,1)])
set(gca,'FontSize',20);

% 

first_landmark_position = 32; % This is set to be the position (in bins) where the motor first starts to move. Whisker touches don't happen until at least bin 32.
second_landmark_position = 72; % This is set to be the position (in bins) where the motor first starts to move. Whisker touches don't happen until at least bin 72.

[amap_lcs_sorted,sort_ind] = sb.sort.rasterMapByPeak(amap_lcs_on,'normalize_type','zscore');
[~,first_peaks] = max(amap_lcs_sorted');

% Sort session and indices based on the sorted map

first_peak_cells = first_peaks<first_landmark_position & first_peaks>16;
first_peak_cells = find(first_peak_cells);

% Find first cell where we are within the second landmark regime
first_cell_second_landmark = find(first_peaks>50);
first_cell_second_landmark = first_cell_second_landmark(1);

[~,second_peaks] = max(amap_lcs_sorted(first_cell_second_landmark:end,:)');

second_peak_cells = second_peaks<second_landmark_position;
second_peak_cells = find(second_peak_cells)+first_cell_second_landmark;
%sb.plot.rasterMaps(amap_lcs_sorted([first_peak_cells,second_peak_cells],:))

% Save rmaps for 50 laps of all predictive cells 
% predictive_rmaps = zeros(50,105);
predictive_rmaps_introduction= [];
for c = 1:length(first_peak_cells)
   predictive_rmaps_introduction = [predictive_rmaps_introduction; amap_lcs_sorted(first_peak_cells(c),:)];
end

for c = 1:length(second_peak_cells)
   predictive_rmaps_introduction = [predictive_rmaps_introduction; amap_lcs_sorted(second_peak_cells(c),:)];
end

predictive_rmaps = predictive_rmaps(:,:,2:end);

predictive_cells_session = session(first_peak_cells);
predictive_cells_inds = indices(first_peak_cells);
predictive_cells_session = [predictive_cells_session, session(second_peak_cells)];
predictive_cells_inds = [predictive_cells_inds, indices(second_peak_cells)];

for i = 1:length(mData_introduction)
    predictive_cells_introduction(i) = 100*length(find(predictive_cells_session == i))/length(unique([mData_introduction(i).rmaps_on.lcs.ind]));
end

nanmean(predictive_cells_introduction)
nanstd(predictive_cells_introduction)/sqrt(length(mData_introduction))



predictive_rmaps_introduction_prct = 100*size(predictive_rmaps_introduction,1)./size(amap_lcs_sorted,1);

line_width = 1;
line_color = 'r';

% Top
line([0,105],[first_peak_cells(1),first_peak_cells(1)],'Color',line_color,'LineWidth',line_width);
line([0,105],[first_peak_cells(end),first_peak_cells(end)],'Color',line_color,'LineWidth',line_width);
line([1,1],[first_peak_cells(1),first_peak_cells(end)],'Color',line_color,'LineWidth',line_width);
line([105,105],[first_peak_cells(1),first_peak_cells(end)],'Color',line_color,'LineWidth',line_width);

% Bottom
line([0,105],[second_peak_cells(1),second_peak_cells(1)],'Color',line_color,'LineWidth',line_width);
line([0,105],[second_peak_cells(end),second_peak_cells(end)],'Color',line_color,'LineWidth',line_width);
line([1,1],[second_peak_cells(1),second_peak_cells(end)],'Color',line_color,'LineWidth',line_width);
line([105,105],[second_peak_cells(1),second_peak_cells(end)],'Color',line_color,'LineWidth',line_width);

% Landmark onset
line([32,32],[0,2000],'Color','w','LineWidth',line_width)
line([72,72],[0,2000],'Color','w','LineWidth',line_width)
xticks([])

% Set correct size of the figure
set(gca,'FontSize',16)
set(gcf,'renderer', 'painters', 'Position', [2200,100,400,700])



