% Here the prediction is done using 0.5 s long windows but I do the
% prediction on every time bin that I record. So basically 0.5s before the
% first bin, then 0.5s before the second bin etc .. this is in difference
% to the way mao does it which is 0.5s bins but non overlapping. 


%% ----- Load data -----
% Clean up workspace
clearvars -except mData sData data workspace
close all;

% Parameters
rng(123) % set random number generator seed
binning = 6; % binning is determined by finding the bin size that gives the minimum mean decoding error across all sessions

% Init
decoding = [];

for session = [1:10]
    
    sData = mData(session).sData;
    
    % Include only place cels
    cells_to_include = [mData(session).rmaps.pcs.ind];
    
    if session == 1
        decoding = decodeSession(sData,cells_to_include);
    else
        decoding(session) = decodeSession(sData,cells_to_include);
    end
end

% disp('done')
%
params.first_landmark_bins = [23:43];
params.second_landmark_bins = 63:83;
otherPos = setdiff(1:105, [params.first_landmark_bins,params.second_landmark_bins]);
% Plot mean decoding error across all sessions
mean_decoding_errors_ON = [];
mean_decoding_errors_OFF = [];
lap_by_lap_decoding_errors_on = [];
lap_by_lap_decoding_errors_off = [];

for s = 1:length(decoding)
    
    % Mean for the whole session
    mean_decoding_errors_ON = [mean_decoding_errors_ON,nanmean(decoding(s).decoding_error_ON)];
    mean_decoding_errors_OFF = [mean_decoding_errors_OFF,nanmean(decoding(s).decoding_error_OFF)];
    
    % Lap by lap erors
    lap_by_lap_decoding_errors_on = [lap_by_lap_decoding_errors_on; decoding(s).decoding_error_at_landmark_on];
    lap_by_lap_decoding_errors_off = [lap_by_lap_decoding_errors_off; decoding(s).decoding_error_at_landmark_off];
    
    idxFlap1= ismember(decoding(s).pos_val_binned, params.first_landmark_bins);
    idxFlap2= ismember(decoding(s).pos_val_binned, params.second_landmark_bins);
    idxOther= ismember(decoding(s).pos_val_binned, otherPos);

    mean_InFlap1(s)    = nanmean(decoding(s).decoding_error_ON(idxFlap1)); 
    mean_InFlap2(s)    = nanmean(decoding(s).decoding_error_ON(idxFlap2));
    mean_Other(s)      = nanmean(decoding(s).decoding_error_ON(idxOther));
    
    mean_BothFlaps(s)    = nanmean([decoding(s).decoding_error_ON(idxFlap1),decoding(s).decoding_error_ON(idxFlap2)]); 

end


for s = 1:length(decoding)
    confMatrix_OFF(s,:,:) = classificationAccuracy(decoding(s).predicted_pos_test, decoding(s).pos_test_binned,unique(decoding(s).pos_test_binned));
    confMatrix_ON(s,:,:) = classificationAccuracy(decoding(s).predicted_pos_val, decoding(s).pos_val_binned,unique(decoding(s).pos_val_binned));
end

decoding_Obj = [];

for session = [1:10]
    
    sData = mData(session).sData;
    
    % Include only place cels
    cells_to_include = [mData(session).rmaps.lcs.ind,mData(session).rmaps.tcs.ind];
    
    if session == 1
        decoding_Obj = decodeSession(sData,cells_to_include);
    else
        decoding_Obj(session) = decodeSession(sData,cells_to_include);
    end
end


% Plot mean decoding error across all sessions
mean_decoding_errors_ON_Obj = [];
mean_decoding_errors_OFF_Obj = [];
lap_by_lap_decoding_errors_on_Obj = [];
lap_by_lap_decoding_errors_off_Obj = [];



for s = 1:length(decoding_Obj)
    
    % Mean for the whole session
    mean_decoding_errors_ON_Obj = [mean_decoding_errors_ON_Obj,nanmean(decoding_Obj(s).decoding_error_ON)];
    mean_decoding_errors_OFF_Obj = [mean_decoding_errors_OFF_Obj,nanmean(decoding_Obj(s).decoding_error_OFF)];
    
    % Lap by lap erors
    lap_by_lap_decoding_errors_on_Obj = [lap_by_lap_decoding_errors_on_Obj; decoding_Obj(s).decoding_error_at_landmark_on];
    lap_by_lap_decoding_errors_off_Obj = [lap_by_lap_decoding_errors_off_Obj; decoding_Obj(s).decoding_error_at_landmark_off];
    
    idxFlap1= ismember(decoding_Obj(s).pos_val_binned, params.first_landmark_bins);
    idxFlap2= ismember(decoding_Obj(s).pos_val_binned, params.second_landmark_bins);
    idxOther= ismember(decoding_Obj(s).pos_val_binned, otherPos);

    mean_InFlap1_Obj(s)    = nanmean(decoding_Obj(s).decoding_error_ON(idxFlap1)); 
    mean_InFlap2_Obj(s)    = nanmean(decoding_Obj(s).decoding_error_ON(idxFlap2)); 
    mean_Other_Obj(s)      = nanmean(decoding_Obj(s).decoding_error_ON(idxOther));
    mean_BothFlaps_Obj(s)    = nanmean([decoding_Obj(s).decoding_error_ON(idxFlap1),decoding_Obj(s).decoding_error_ON(idxFlap2)]); 


end


for s = 1:length(decoding_Obj)
    confMatrix_OFF_Obj(s,:,:) = classificationAccuracy(decoding_Obj(s).predicted_pos_test, decoding_Obj(s).pos_test_binned,unique(decoding_Obj(s).pos_test_binned));
    confMatrix_ON_Obj(s,:,:) = classificationAccuracy(decoding_Obj(s).predicted_pos_val, decoding_Obj(s).pos_val_binned,unique(decoding_Obj(s).pos_val_binned));
end
       

%% Show decoding error as a funciton of position with error shading the SEM
figure(5); clf;

% Landmark ON
x = 1:105; % bins
y = sb.generate.averagedTuningFromRasterMaps(lap_by_lap_decoding_errors_on);
y = smoothdata(y,'movmean',4);
SEM = nanstd(lap_by_lap_decoding_errors_on)/sqrt(size(lap_by_lap_decoding_errors_on,1));
s = shadedErrorBar(x,y,SEM);
s.mainLine.LineWidth = 2;

hold on;

% Landmark OFF
x = 1:105; % bins
y = sb.generate.averagedTuningFromRasterMaps(lap_by_lap_decoding_errors_off);
y = smoothdata(y,'movmean',4);
SEM = nanstd(lap_by_lap_decoding_errors_off)/sqrt(size(lap_by_lap_decoding_errors_off,1));
s = shadedErrorBar(x,y,SEM);
s.patch.FaceColor = [0,0,1];
s.mainLine.Color = [0,0,1];
s.mainLine.LineWidth = 2;
s.patch.FaceAlpha = 0.3;

xlim([1,105]);
xticks([1,105]);
xticklabels({0,157})
xlabel('Position (cm)');

yticks([0,20]);
ylabel('Decoding error (cm)');
ylim([0,24])
set(gca,'FontSize',16);

% Draw landmark lines
line([32,32],[0,100],'Color','r');
line([72,72],[0,100],'Color','r');
line([92,92],[0,100],'Color','k');

legend({'Landmark present','Landmark omitted'})


% Set figure size
set(gcf,'renderer', 'painters', 'Position', [200,1000,1200,500])



mean_confMatrix_ON  = nanmean(confMatrix_ON(:,:,:),1);
mean_confMatrix_OFF  = nanmean(confMatrix_OFF(:,:,:),1);

mean_confMatrix_ON_Obj  = nanmean(confMatrix_ON_Obj(:,:,:),1);
mean_confMatrix_OFF_Obj   = nanmean(confMatrix_OFF_Obj(:,:,:),1);

figure()
subplot(1,3,1)
imagesc(flipud(reshape(mean_confMatrix_ON,105,105)))
ax = gca;
colorbar
ax.XColor = [0 0 0];
ax.YColor = [0 0 0];
ax.XLabel.Color = [0 0 0];
ax.YLabel.Color = [0 0 0];
set(gca,'xTick',[1 33.33 66.66 105], 'xTickLabel',{'1','50','100','157'})
set(gca,'yTick',[1 105-66.66 105-33.33  105], 'yTickLabel',{'157','100','50','1'});
set(gca,'FontSize', 18, 'FontName','Gotham');
colorbar()
xline(33,'Color','w','LineWidth',1)
xline(73,'Color','w','LineWidth',1)
yline(105-33,'Color','w','LineWidth',1)
yline(105-73,'Color','w','LineWidth',1)
xlabel('Position (cm)')
ylabel('Predicted Position (cm)')
title('Position tuned cells')

subplot(1,3,2)
imagesc(flipud(reshape(mean_confMatrix_ON_Obj,105,105)))
ax = gca;
colorbar
ax.XColor = [0 0 0];
ax.YColor = [0 0 0];
ax.XLabel.Color = [0 0 0];
ax.YLabel.Color = [0 0 0];
set(gca,'xTick',[1 33.33 66.66 105], 'xTickLabel',{'1','50','100','157'})
set(gca,'yTick',[1 105-66.66 105-33.33  105], 'yTickLabel',{'157','100','50','1'});
set(gca,'FontSize', 18, 'FontName','Gotham');
colorbar()
xline(33,'Color','w','LineWidth',1)
xline(73,'Color','w','LineWidth',1)
yline(105-33,'Color','w','LineWidth',1)
yline(105-73,'Color','w','LineWidth',1)
xlabel('Position (cm)')
ylabel('Predicted Position (cm)')
title('Object location cells')



subplot(1,3,3)
bar([nanmean(mean_decoding_errors_ON),nanmean(mean_decoding_errors_ON_Obj)],'FaceColor',[0.5 0.5 0.5])
hold on 
errorbar(1:2,[nanmean(mean_decoding_errors_ON),nanmean(mean_decoding_errors_ON_Obj)],...
    [nanstd(mean_decoding_errors_ON)/sqrt(length(mean_decoding_errors_ON_Obj)),nanstd(mean_decoding_errors_ON_Obj)/sqrt(length(mean_decoding_errors_ON_Obj))],'k','LineStyle','none',...
    'LineWidth',1)
ax.XColor = [0 0 0];
ax.YColor = [0 0 0];
ax.XLabel.Color = [0 0 0];
ax.YLabel.Color = [0 0 0];
xticks([1  2]); xticklabels({'Position tuned cells','Object location cells'})
set(gca,'FontSize', 18, 'FontName','Gotham');
yticks([0:10:40])
ylabel({'Decoding error (cm)'})
yline(39.25)
xtickangle(45)


% pval = ranksum(mean_decoding_errors_ON_Obj,mean_decoding_errors_ON);
ranksum(mean_decoding_errors_ON,ones(10,1)*39.25,'tail','left')
ranksum(mean_decoding_errors_ON_Obj,ones(10,1)*39.25,'tail','left')

%
% figure()
% subplot(2,1,1)
% y = nanmean(mean_confMatrix_ON(:,:,params.first_landmark_bins),3);
% SEM = nanstd(reshape(mean_confMatrix_ON(:,:,params.first_landmark_bins),105,length(params.first_landmark_bins))')./sqrt(length(params.first_landmark_bins));
% y = smoothdata(y,'movmean',4);
% s = shadedErrorBar(1:105,y,SEM);
% s.patch.FaceColor = [0,0,0];
% s.mainLine.Color = [0,0,0];
% s.mainLine.LineWidth = 2;
% s.patch.FaceAlpha = 0.3;
% hold on
% 
% y = nanmean(mean_confMatrix_ON_Obj(:,:,params.first_landmark_bins),3);
% SEM = nanstd(reshape(mean_confMatrix_ON_Obj(:,:,params.first_landmark_bins),105,length(params.first_landmark_bins))')./sqrt(length(params.first_landmark_bins));
% y = smoothdata(y,'movmean',4);
% s = shadedErrorBar(1:105,y,SEM);
% s.patch.FaceColor = [0,0,1];
% s.mainLine.Color = [0,0,1];
% s.mainLine.LineWidth = 2;
% s.patch.FaceAlpha = 0.3;
% xlim([0 105])
% yticks([0  0.02 0.04])
% yticklabels(100*[0  0.02 0.04])
% ylabel('Prediction position 1st flap(%)')
% 
% xlim([0 105])
% yticks([0  0.02 0.04 0.06 0.08])
% yticklabels(100*[0  0.02 0.04 0.06 0.08])
% ylabel('Prediction position 1st flap(%)')
% 
% ax = gca;
% ax.XColor = [0 0 0];
% ax.YColor = [0 0 0];
% ax.XLabel.Color = [0 0 0];
% ax.YLabel.Color = [0 0 0];
% 
% xtix = get(gca, 'XTick');
% set(gca, 'XTick',xtix, 'XTickLabel',xtix*1.5)
% xticks([0,33,73,105])
% xticklabels({'0','50','110','157'})
% xlabel('Position (cm)');
% %title('Stimulus ON');
% set(gca,'FontSize',18);
% xline(33,'Color','k','LineWidth',1)
% xline(73,'Color','k','LineWidth',1)
% 
% subplot(2,1,2)
% y = nanmean(mean_confMatrix_ON(:,:,params.second_landmark_bins),3);
% SEM = nanstd(reshape(mean_confMatrix_ON(:,:,params.second_landmark_bins),105,length(params.second_landmark_bins))')./sqrt(length(params.second_landmark_bins));
% y = smoothdata(y,'movmean',4);
% s = shadedErrorBar(1:105,y,SEM);
% s.mainLine.LineWidth = 2;
% s.patch.FaceAlpha = 0.3;
% hold on
% 
% y = nanmean(mean_confMatrix_ON_Obj(:,:,params.second_landmark_bins),3);
% SEM = nanstd(reshape(mean_confMatrix_ON_Obj(:,:,params.second_landmark_bins),105,length(params.second_landmark_bins))')./sqrt(length(params.second_landmark_bins));
% y = smoothdata(y,'movmean',4);
% s = shadedErrorBar(1:105,y,SEM);
% s.patch.FaceColor = [0,0,1];
% s.mainLine.Color = [0,0,1];
% s.mainLine.LineWidth = 2;
% s.patch.FaceAlpha = 0.3;
% 
% xlim([0 105])
% yticks([0  0.02 0.04 0.06 0.08])
% yticklabels(100*[0  0.02 0.04 0.06 0.08])
% ylabel('Prediction position 2nd flap(%)')
% 
% ax = gca;
% ax.XColor = [0 0 0];
% ax.YColor = [0 0 0];
% ax.XLabel.Color = [0 0 0];
% ax.YLabel.Color = [0 0 0];
% 
% xtix = get(gca, 'XTick');
% set(gca, 'XTick',xtix, 'XTickLabel',xtix*1.5)
% xticks([0,33,73,105])
% xticklabels({'0','50','110','157'})
% xlabel('Position (cm)');
% %title('Stimulus ON');
% set(gca,'FontSize',18);
% xline(33,'Color','k','LineWidth',1)
% xline(73,'Color','k','LineWidth',1)




% prediction accuracy first landmark vs 2nd landmark

% figure()
% 
% bar([1 3 5],[nanmean(mean_InFlap1),nanmean(mean_InFlap2),nanmean(mean_Other)],'FaceColor',[0.5 0.5 0.5],'BarWidth',0.3)
% hold on
% errorbar([1 3 5],[nanmean(mean_InFlap1),nanmean(mean_InFlap2),nanmean(mean_Other)],...
%     [nanstd(mean_InFlap1)/sqrt(length(mean_InFlap1)),nanstd(mean_InFlap2)/sqrt(length(mean_InFlap1)),nanstd(mean_Other)/sqrt(length(mean_Other))],'k','LineStyle','none','LineWidth',1)
% hold on
% bar([2 4 6],[nanmean(mean_InFlap1_Obj),nanmean(mean_InFlap2_Obj),nanmean(mean_Other_Obj)],'FaceColor','b','BarWidth',0.3)
% hold on
% errorbar([2 4  6],[nanmean(mean_InFlap1_Obj),nanmean(mean_InFlap2_Obj),nanmean(mean_Other_Obj)],...
%     [nanstd(mean_InFlap1_Obj)/sqrt(length(mean_InFlap1_Obj)),nanstd(mean_InFlap2_Obj)/sqrt(length(mean_InFlap1_Obj)),nanmean(mean_Other_Obj)/sqrt(length(mean_Other_Obj))],'k','LineStyle','none','LineWidth',1)
% ax = gca;
% ax.XColor = [0 0 0];
% ax.YColor = [0 0 0];
% ax.XLabel.Color = [0 0 0];
% ax.YLabel.Color = [0 0 0];
% xticks([1.5  3.5  5.5]); xticklabels({'Other Positions','1st flap','2nd flap'})
% set(gca,'FontSize', 18, 'FontName','Gotham');
% ylabel({'Decoding error (cm)'})
% xtickangle(45)
% xlim([0.3 6.8])
% 
% pval(1) = signrank(mean_Other,mean_Other_Obj,'tail','left');
% pval(2) = signrank(mean_InFlap1,mean_InFlap1_Obj,'tail','left');
% pval(3) = signrank(mean_InFlap2,mean_InFlap2_Obj,'tail','left');
% pval(4) = signrank(mean_Other_Obj,mean_InFlap1_Obj,'tail','left');
% pval(5) = signrank(mean_Other_Obj,mean_InFlap2_Obj,'tail','left');
% 
% box off
% groups={[1,2],[3,4],[5,6],[2 4],[2 6]};
% sigstar(groups,pval);
% 
% yline(39.25)

% 
% figure()
% 
% bar(1,nanmean(mean_BothFlaps),'FaceColor',[0.5 0.5 0.5],'BarWidth',0.5); hold on
% errorbar(1,nanmean(mean_BothFlaps),nanstd(mean_BothFlaps)./sqrt(length(mean_BothFlaps)),'k','LineStyle','none','LineWidth',1)
% bar(2,nanmean(mean_BothFlaps_Obj),'FaceColor','b','BarWidth',0.5)
% errorbar(2,nanmean(mean_BothFlaps_Obj),nanstd(mean_BothFlaps_Obj)./sqrt(length(mean_BothFlaps)),'k','LineStyle','none','LineWidth',1)
% pval = [];
% pval(1) = signrank(mean_BothFlaps,mean_BothFlaps_Obj,'tail','left');
% groups={[1,2]};
% sigstar(groups,pval);
% 
% ax.XColor = [0 0 0];
% ax.YColor = [0 0 0];
% ax.XLabel.Color = [0 0 0];
% ax.YLabel.Color = [0 0 0];
% xticks([1 2]); xticklabels({'Position tuned cells','Object location cells'})
% xtickangle(45)
% ylabel({'Decoding error of object positions(cm)'})
% 
% yline(39.25)
%% ---- HELPER FUNCTION ----- 
    function decoding_results = decodeSession(sData,rois_to_use)

    if nargin < 2
        rois_to_use = [];
    end
    
    % Parameters
    n_bins_to_pool = 6; % How many time bins should be binned for the Bayesian model. This was determined by finding the bin size that minimized the mean decoding error across all sessions. 6 bins is ca 0.2 seconds.
    tau = n_bins_to_pool/31; % 31 is the fps of the imaging. Tau is used for the Bayesian model
    rng(123) % set random number generator seed
    n_spatial_bins = 105; % number of bins along the linear track
    
    % Init vars
    decoding_results = struct();
    
    % Extract data
    deconv = double(sData.imdata.roiSignals(2).deconv);
    pos = sData.behavior.wheelPosDsBinned;
    lap = sData.behavior.wheelLapDsBinned;
    speed = sData.behavior.runSpeedDs;

    % Exclude silent cells
%     deconv = deconv(find(sum(deconv,2)),:);

    % Remove all bins where running speed is < 1 cm.
    bins_to_include = find([speed<1] == 0);
    deconv = deconv(:,bins_to_include);
    pos = pos(bins_to_include);
    lap = lap(bins_to_include);

    % Add non-zero noise to deconv
    deconv = deconv + randi(10,size(deconv,1),size(deconv,2))*0.0000000000000001;

    % Smooth the deconv signal (Mao does this in his code, but don't write so in the paper) 
    deconv = smoothdata(deconv,2,'gaussian',5);
    
    % Include only a subset of cells
    if ~isempty(rois_to_use)
       
        % Subset the deconv matrix
        deconv = deconv(rois_to_use,:);
        
    end

    % Shift lap to start on 1, not -1 or 0.
    lap = lap + (abs(min(lap))*2);

    % Find laps where the landmark motor is present (on) or omitted (off)
    laps_landmark_on  = find(sData.landmarks.motor.on(:,1) == 1);
    laps_landmark_off = find(sData.landmarks.motor.on(:,1) == 0);

    % Split into train and test set (odd vs even laps)
    laps_landmark_on_even = laps_landmark_on(find(mod(laps_landmark_on,2) == 0));
    laps_landmark_on_odd  = laps_landmark_on(find(mod(laps_landmark_on,2) == 1));

    % Find indices for the OFF laps, ON (even) and ON (odd) laps
    off_laps     = find(ismember(lap,laps_landmark_off));
    on_odd_laps  = find(ismember(lap,laps_landmark_on_odd));
    on_even_laps = find(ismember(lap,laps_landmark_on_even));

    % Find the shortest number of bins between the even laps and off laps
    % to balance the validation and test set
    max_number_of_inds = min([length(on_even_laps),length(off_laps)]);
    on_even_laps = on_even_laps(1:max_number_of_inds);
    off_laps = off_laps(1:max_number_of_inds);
    
    % Divide the data into train, val, test
    deconv_train = deconv(:,on_odd_laps);
    deconv_val  = deconv(:,on_even_laps);
    deconv_test = deconv(:,off_laps);
    pos_train = pos(on_odd_laps);
    pos_test  = pos(off_laps);
    pos_val   = pos(on_even_laps);
    lap_train = lap(on_odd_laps);
    lap_test  = lap(off_laps);
    lap_val   = lap(on_even_laps);
    
    % ----- Compute position tuning curves (1D) -----
    % Compute trial averaged and occupancy-normalized position tuning curves
    % for all neurons, using the training dataset

    % Occupancy 
    occupancy = zeros(1, n_spatial_bins);
    for i = 1:length(pos_train)
        
       % Add one to each position in the occupancy array
       occupancy(pos_train(i)) = occupancy(pos_train(i)) + 1;
       
    end

    position_tuning_curves = [];

    % Position tuning curves
    % For each cell
    for c = 1:size(deconv_train,1)
        sig = deconv_train(c,:);

        % Tuning map of the current cell
        cell_tuning = zeros(1,n_spatial_bins);

        % For each datapoint in the signal add the signal to corresponding
        % position
        for i = 1:length(sig)
            cell_tuning(pos_train(i)) = cell_tuning(pos_train(i)) + sig(i);
        end

        % Normalize the cell tuning with occupancy
        cell_tuning = cell_tuning ./ occupancy;

        % Smooth the cell tuning
        cell_tuning = smoothdata(cell_tuning,'gaussian',3);

        % Add cell tuning to position tuning curve matrix
        if c == 1
            position_tuning_curves = cell_tuning;
        else
            position_tuning_curves = [position_tuning_curves;cell_tuning];
        end
    end


    % Compute the multiplying factor 
    multiply_factor = exp(-tau * sum(position_tuning_curves,1)); % I use dim 1 to sum across all cells

    % ----- Predict position -----

    % Bin the test data of deconv
    deconv_test_binned = [];
    deconv_val_binned = [];
    pos_test_binned = [];
    pos_val_binned = [];

    % Test data (OFF laps)
    i = 1;
    for s = 1:n_bins_to_pool:length(deconv_test)-n_bins_to_pool
        deconv_test_binned(:,i) = nanmean(deconv_test(:,s:s+n_bins_to_pool-1),2);
        pos_test_binned(:,i) = pos_test(:,s+n_bins_to_pool-1);
        lap_test_binned(:,i) = lap_test(:,s+n_bins_to_pool-1);
        i = i+1;
    end

    % Validation data (ON even laps)
    i = 1;
    for s = 1:n_bins_to_pool:length(deconv_val)-n_bins_to_pool
        deconv_val_binned(:,i) = nanmean(deconv_val(:,s:s+n_bins_to_pool-1),2);
        pos_val_binned(:,i) = pos_val(:,s+n_bins_to_pool-1);
        lap_val_binned(:,i) = lap_val(:,s+n_bins_to_pool-1);
        i = i+1;
    end

    % Predict position of test data for each time bin
    predicted_pos_test = [];
    for t = 1:size(deconv_test_binned,2)

        f_pos = position_tuning_curves(:,:);
        mean_activity = deconv_test_binned(:,t);
        result = prod(f_pos.^mean_activity,1);
        result = result / sum(result(:));
        result = result.* multiply_factor;

        [~,predicted_pos_test(t)] = max(result);

    end

    % Predict position of validation data for each time bin
    predicted_pos_val = [];
    for t = 1:size(deconv_val_binned,2)

        f_pos = position_tuning_curves(:,:);
        mean_activity = deconv_val_binned(:,t);
        result = prod(f_pos.^mean_activity,1);
        result = result / sum(result(:));
        result = result.* multiply_factor;

        [~,predicted_pos_val(t)] = max(result);
    end
    

    % Compute the decoding error
    decoding_error_test = abs(pos_test_binned-predicted_pos_test);
    decoding_error_val  = abs(pos_val_binned-predicted_pos_val);

    % Correct for circularity
    need_correcting = find(decoding_error_test>(n_spatial_bins/2));
    decoding_error_test(need_correcting) = n_spatial_bins-decoding_error_test(need_correcting);
    mean_decoding_error_test = nanmean(decoding_error_test)*1.5;

    need_correcting = find(decoding_error_val>(n_spatial_bins/2));
    decoding_error_val(need_correcting) = n_spatial_bins-decoding_error_val(need_correcting);
    mean_decoding_error_val = nanmean(decoding_error_val)*1.5;

    % For each lap get the decoding error average between position 50 and 60
    % cm.
    
    % Compute the average error after landmark when not omitted
    decoding_error_at_landmark_on = [];
    lap_num = 1;
    for lap = unique(lap_val_binned)
        
        % Compute for each bin
        for current_bin = 1:105
            decoding_error_at_landmark_on(lap_num,current_bin) = nanmean(decoding_error_val([lap_val_binned == lap] & ismember(pos_val_binned,current_bin)));
        end
        lap_num = lap_num + 1;
    end

    % Compute the average error at after landmark when landmark is omitted
    decoding_error_at_landmark_off = [];
    lap_num = 1;
    for lap = unique(lap_test_binned)
        for current_bin = 1:105
            decoding_error_at_landmark_off(lap_num,current_bin) = nanmean(decoding_error_test([lap_test_binned == lap] & ismember(pos_test_binned,current_bin)));
        end
        lap_num = lap_num + 1;
    end

    % Find the common minimum number of laps to include to make the lengths
    % equal
    max_laps = min([size(decoding_error_at_landmark_on,1),size(decoding_error_at_landmark_off,1)]);
    decoding_error_at_landmark_on = decoding_error_at_landmark_on(1:max_laps,:);
    decoding_error_at_landmark_off = decoding_error_at_landmark_off(1:max_laps,:);
   
    
    % Add data to output
    decoding_results.decoding_error_ON = decoding_error_val*1.5;
    decoding_results.decoding_error_OFF = decoding_error_test*1.5;
   
    decoding_results.pos_val_binned = pos_val_binned;
    decoding_results.predicted_pos_val = predicted_pos_val;
    
    decoding_results.pos_test_binned = pos_test_binned;
    decoding_results.predicted_pos_test = predicted_pos_test;
    
    decoding_results.pos_test_original = pos_test;
    
    decoding_results.pos_train = pos_train;
    decoding_results.pos_val = pos_val;
    decoding_results.pos_test = pos_test;
    
    % Lap by lap error
    decoding_results.decoding_error_at_landmark_on = decoding_error_at_landmark_on;
    decoding_results.decoding_error_at_landmark_off = decoding_error_at_landmark_off;
    
    
    end

    function norm_confMatrix = classificationAccuracy(predClass, actualClass,classes)

        classAcc = zeros(length(classes),1);

        for i = 1: length(classes)
            idx                   = actualClass == classes(i);
            predClasstemp         = predClass(idx);
            corrPred(i)           = numel(find(predClasstemp == classes(i)));
            classAcc(i)           = corrPred(i)/length(find(idx==1));

            norm_confMatrix(i,i)  = corrPred(i)/length(find(idx==1));

            falsePred = zeros(length(classes),1);
            for j = 1: length(classes)
                if j == i
                else
                    falsePred(j) = length(find(predClasstemp == classes(j)));
                    norm_confMatrix(j,i) = falsePred(j)/length(find(idx==1));
                end
            end


        end
    end
    
    
    
    function varargout=sigstar(groups,stats,nosort)
    % sigstar - Add significance stars to bar charts, boxplots, line charts, etc,
    %
    % H = sigstar(groups,stats,nsort)
    %
    % Purpose
    % Add stars and lines highlighting significant differences between pairs of groups. 
    % The user specifies the groups and associated p-values. The function handles much of 
    % the placement and drawing of the highlighting. Stars are drawn according to:
    %   * represents p<=0.05
    %  ** represents p<=1E-2
    % *** represents p<=1E-3
    %
    %
    % Inputs
    % groups - a cell array defining the pairs of groups to compare. Groups defined 
    %          either as pairs of scalars indicating locations along the X axis or as 
    %          strings corresponding to X-tick labels. Groups can be a mixture of both 
    %          definition types.
    % stats -  a vector of p-values the same length as groups. If empty or missing it's 
    %          assumed to be a vector of 0.05s the same length as groups. Nans are treated
    %          as indicating non-significance.
    % nsort -  optional, 0 by default. If 1, then significance markers are plotted in 
    %          the order found in groups. If 0, then they're sorted by the length of the 
    %          bar.
    %
    % Outputs
    % H - optionally return handles for significance highlights. Each row is a different
    %     highlight bar. The first column is the line. The second column is the text (stars).
    %     
    %
    % Examples
    % 1. 
    % bar([5,2,1.5])
    % sigstar({[1,2], [1,3]})
    %
    % 2. 
    % bar([5,2,1.5])
    % sigstar({[2,3],[1,2], [1,3]},[nan,0.05,0.05])
    %
    % 3.  **DOESN'T WORK IN 2014b**
    % R=randn(30,2);
    % R(:,1)=R(:,1)+3;
    % boxplot(R)
    % set(gca,'XTick',1:2,'XTickLabel',{'A','B'})
    % H=sigstar({{'A','B'}},0.01);
    % ylim([-3,6.5])
    % set(H,'color','r')
    %
    % 4. Note the difference in the order with which we define the groups in the 
    %    following two cases. 
    % x=[1,2,3,2,1];
    % subplot(1,2,1)
    % bar(x)
    % sigstar({[1,2], [2,3], [4,5]})
    % subplot(1,2,2)
    % bar(x)
    % sigstar({[2,3],[1,2], [4,5]})
    %
    % ALSO SEE: demo_sigstar
    %
    % KNOWN ISSUES:
    % 1. Algorithm for identifying whether significance bar will overlap with 
    %    existing plot elements may not work in some cases (see line 277)
    % 2. Bars may not look good on exported graphics with small page sizes.
    %    Simply increasing the width and height of the graph with the 
    %    PaperPosition property of the current figure should fix things.
    %
    % Rob Campbell - CSHL 2013
    %Input argument error checking
    %If the user entered just one group pair and forgot to wrap it in a cell array 
    %then we'll go easy on them and wrap it here rather then generate an error
    if ~iscell(groups) & length(groups)==2
        groups={groups};
    end
    if nargin<2 
        stats=repmat(0.05,1,length(groups));
    end
    if isempty(stats)
        stats=repmat(0.05,1,length(groups));
    end
    if nargin<3
        nosort=0;
    end
    %Check the inputs are of the right sort
    if ~iscell(groups)
        error('groups must be a cell array')
    end
    if ~isvector(stats)
        error('stats must be a vector')
    end
    if length(stats)~=length(groups)
        error('groups and stats must be the same length')
    end
    %Each member of the cell array groups may be one of three things:
    %1. A pair of indices.
    %2. A pair of strings (in cell array) referring to X-Tick labels
    %3. A cell array containing one index and one string
    %
    % For our function to run, we will need to convert all of these into pairs of
    % indices. Here we loop through groups and do this. 
    xlocs=nan(length(groups),2); %matrix that will store the indices 
    xtl=get(gca,'XTickLabel');  
    for ii=1:length(groups)
        grp=groups{ii};
        if isnumeric(grp)
            xlocs(ii,:)=grp; %Just store the indices if they're the right format already
        elseif iscell(grp) %Handle string pairs or string/index pairs
            if isstr(grp{1})
                a=strmatch(grp{1},xtl);
            elseif isnumeric(grp{1})
                a=grp{1};
            end
            if isstr(grp{2})
                b=strmatch(grp{2},xtl);
            elseif isnumeric(grp{2})
                b=grp{2};
            end
            xlocs(ii,:)=[a,b];
        end
        %Ensure that the first column is always smaller number than the second
        xlocs(ii,:)=sort(xlocs(ii,:));
    end
    %If there are any NaNs we have messed up. 
    if any(isnan(xlocs(:)))
        error('Some groups were not found')
    end
    %Optionally sort sig bars from shortest to longest so we plot the shorter ones first
    %in the loop below. Usually this will result in the neatest plot. If we waned to 
    %optimise the order the sig bars are plotted to produce the neatest plot, then this 
    %is where we'd do it. Not really worth the effort, though, as few plots are complicated
    %enough to need this and the user can define the order very easily at the command line. 
    if ~nosort
        [~,ind]=sort(xlocs(:,2)-xlocs(:,1),'ascend');
        xlocs=xlocs(ind,:);groups=groups(ind);
        stats=stats(ind);
    end
    %-----------------------------------------------------
    %Add the sig bar lines and asterisks 
    holdstate=ishold;
    hold on
    H=ones(length(groups),2); %The handles will be stored here
    y=ylim;
    yd=myRange(y)*0.05; %separate sig bars vertically by 5% 
    for ii=1:length(groups)
        thisY=findMinY(xlocs(ii,:))+yd;
        H(ii,:)=makeSignificanceBar(xlocs(ii,:),thisY,stats(ii));
    end
    %-----------------------------------------------------
    %Now we can add the little downward ticks on the ends of each line. We are
    %being extra cautious and leaving this it to the end just in case the y limits
    %of the graph have changed as we add the highlights. The ticks are set as a
    %proportion of the y axis range and we want them all to be the same the same
    %for all bars.
    yd=myRange(ylim)*0.03; %Ticks are 1% of the y axis range
    for ii=1:length(groups)
        y=get(H(ii,1),'YData');
        y(1)=y(1)-yd;
        y(4)=y(4)-yd;   
        set(H(ii,1),'YData',y)
    end
    %Be neat and return hold state to whatever it was before we started
    if ~holdstate
        hold off
    elseif holdstate
        hold on
    end
    %Optionally return the handles to the plotted significance bars (first column of H)
    %and asterisks (second column of H).
    if nargout>0
        varargout{1}=H;
    end
end %close sigstar

function H=makeSignificanceBar(x,y,p)
    %makeSignificanceBar produces the bar and defines how many asterisks we get for a 
    %given p-value
    if p<=1E-3
        stars='***'; 
    elseif p<=1E-2
        stars='**';
    elseif p<=0.05
        stars='*';
    elseif isnan(p)
        stars='n.s.';
    else
        stars='n.s.';
    end
            
    x=repmat(x,2,1);
    y=repmat(y,4,1);
    H(1)=plot(x(:),y,'-k','LineWidth',1.5,'Tag','sigstar_bar');
    %Increase offset between line and text if we will print "n.s."
    %instead of a star. 
%     if ~isnan(p)
%         offset=0.005;
%     else
        offset=0.02;
%     end
    starY=mean(y)+myRange(ylim)*offset;
    H(2)=text(mean(x(:)),double(starY),stars,...
        'HorizontalAlignment','Center',...
        'BackGroundColor','none',...
        'Tag','sigstar_stars','FontSize',16);
    Y=ylim;
    if Y(2)<starY
        ylim([Y(1),starY+myRange(Y)*0.05])
    end
end %close makeSignificanceBar

function Y=findMinY(x)
    % The significance bar needs to be plotted a reasonable distance above all the data points
    % found over a particular range of X values. So we need to find these data and calculat the 
    % the minimum y value needed to clear all the plotted data present over this given range of 
    % x values. 
    %
    % This version of the function is a fix from Evan Remington
    oldXLim = get(gca,'XLim');
    oldYLim = get(gca,'YLim');
    axis(gca,'tight')
    
    %increase range of x values by 0.1 to ensure correct y max is used
    x(1)=x(1)-0.1;
    x(2)=x(2)+0.1;
    
    set(gca,'xlim',x) %Matlab automatically re-tightens y-axis
    yLim = get(gca,'YLim'); %Now have max y value of all elements within range.
    Y = max(yLim);
    axis(gca,'normal')
    set(gca,'XLim',oldXLim,'YLim',oldYLim)
end %close findMinY

function rng=myRange(x)
    %replacement for stats toolbox range function
    rng = max(x) - min(x);
end %close myRange

