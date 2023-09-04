% identify object cells for two different types of objects
% object A (at cue_1 positions) is sandpaper
% object B (at cue_2 positions) is felt pad

cue_1 = [42 47 122 128];
cue_1(1) = ceil((cue_1(1)-8)./1.5);
cue_1(3) = ceil((cue_1(3)-8)./1.5);
cue_1(2) = ceil((cue_1(2)+8)./1.5);
cue_1(4) = ceil((cue_1(4)+8)./1.5);

paramsA.first_landmark_bins = cue_1(1):cue_1(2);
paramsA.second_landmark_bins = cue_1(3):cue_1(4);
    

cue_2 = [63 68 100 106];
cue_2(1) = ceil((cue_2(1)-8)./1.5);
cue_2(3) = ceil((cue_2(3)-8)./1.5);
cue_2(2) = ceil((cue_2(2)+8)./1.5);
cue_2(4) = ceil((cue_2(4)+8)./1.5);

paramsB.first_landmark_bins = cue_2(1):cue_2(2);
paramsB.second_landmark_bins = cue_2(3):cue_2(4);

lcsB_rmap = cell(length(data.sessionIDs));
lcsA_rmap = cell(length(data.sessionIDs));
lcsAB_rmap= cell(length(data.sessionIDs));

lcsB_rmap_mean = cell(length(data.sessionIDs));
lcsA_rmap_mean = cell(length(data.sessionIDs));
lcsAB_rmap_mean= cell(length(data.sessionIDs));

lcsA_prct  = NaN(1,length(data.sessionIDs));
lcsB_prct  = NaN(1,length(data.sessionIDs));
lcsAB_prct = NaN(1,length(data.sessionIDs));


for s = 1:length(data.sessionIDs)
    
    lcsA= sb.classify.landmarkCellsByShuffling(mData(s).rmaps.deconv,paramsA);
    lcsB= sb.classify.landmarkCellsByShuffling(mData(s).rmaps.deconv,paramsB);
    
%     if ~exist(fullfile(savedrct,data.sessionIDs{s}));mkdir(fullfile(savedrct,data.sessionIDs{s}));end
%     save(fullfile(savedrct,data.sessionIDs{s},'lcsA.mat'),'lcsA')
%     save(fullfile(savedrct,data.sessionIDs{s},'lcsB.mat'),'lcsB')
    
    mData(s).lcs_sandpaper =  lcsA;
    mData(s).lcs_feltpad   =  lcsB;
    lcsAB = intersect(lcsA,lcsB);
    lcsA  = setdiff(lcsA,lcsAB);
    lcsB  = setdiff(lcsB,lcsAB);
    
    lcsA_rmap{s}        = mData(s).rmaps.deconv(:,:,lcsA);
    lcsA_rmap_mean{s}   = reshape(nanmean(mData(s).rmaps.deconv(:,:,lcsA),1),105,length(lcsA));
    lcsA_prct(s)        = length(lcsA)/size(mData(s).rmaps.deconv,3);
    
    lcsB_rmap_mean{s}   = reshape(nanmean(mData(s).rmaps.deconv(:,:,lcsB),1),105,length(lcsB));
    lcsB_prct(s)        =  length(lcsB)/size(mData(s).rmaps.deconv,3);
    lcsB_rmap{s}        = mData(s).rmaps.deconv(:,:,lcsB);

    
    lcsAB_rmap_mean{s}  = reshape(nanmean(mData(s).rmaps.deconv(:,:,lcsAB),1),105,length(lcsAB));
    lcsAB_prct(s)       = numel(lcsAB)/size(mData(s).rmaps.deconv,3);
    lcsAB_rmap{s}       = mData(s).rmaps.deconv(:,:,lcsAB);

end


% the exact percentages might vary a bit due to the randomneess of the
% shuffling
figFrc = figure();
set(figFrc, 'Units', 'centimeters');
set(figFrc, 'Position', [0 0 18 15]);
scatter(2*ones(length(lcsA_prct),1),100*(lcsA_prct'+lcsAB_prct'),90,'k','LineWidth',1)
hold on
scatter(1*ones(length(lcsB_prct),1),100*(lcsB_prct'+lcsAB_prct'),90,'k','LineWidth',1)
% scatter(3*ones(length(feltPadsandpaperCell_prct),1),100*feltPadsandpaperCell_prct,90,'k','LineWidth',1)
xlim([0.5 2.5])
xticks([1:2])
xticklabels({'felt pad','sandpaper'})
xtickangle(45)
ax = gca;
ax.XColor = [0 0 0];
ax.YColor = [0 0 0];
ax.XLabel.Color = [0 0 0];
ax.YLabel.Color = [0 0 0];
set(gca, 'FontSize',18, 'FontName', 'Arial');
ylabel('Object location cells (%)')

bar([nanmean(100*(lcsB_prct+lcsAB_prct)),nanmean(100*(lcsA_prct+lcsAB_prct))],...
'FaceColor','none')
hold on
errorbar(1:2,[nanmean(100*(lcsB_prct+lcsAB_prct)),nanmean(100*(lcsA_prct+lcsAB_prct))],...
    [nanstd(100*(lcsB_prct+lcsAB_prct))/sqrt(length(lcsB_prct)),nanstd(100*(lcsA_prct+lcsAB_prct))/sqrt(length(lcsA_prct))],...
    'Color','k','LineWidth',1,'LineStyle','none')
ylim([0 35])
box off
xlim([0.5 2.5])
xticks([1:3])
xticklabels({'felt pad','sandpaper'})
xtickangle(45)
ax = gca;
ax.XColor = [0 0 0];
ax.YColor = [0 0 0];
ax.XLabel.Color = [0 0 0];
ax.YLabel.Color = [0 0 0];
set(gca, 'FontSize',18, 'FontName', 'Arial');
ylabel('Object location cells (%)')

saveas(figFrc,fullfile(savedrct,'lanmarkType.fig'),'fig')
saveas(figFrc,fullfile(savedrct,'lanmarkType.png'),'png')
saveas(figFrc,fullfile(savedrct,'lanmarkType.pdf'),'pdf')

figSequence = figure();
set(figSequence, 'Units', 'centimeters');
set(figSequence, 'Position', [0 0 30 15]);
subplot(1,3,1)
lcsBMatrix = cell2mat(lcsB_rmap_mean')';
normedlcsBMatrix    = normalize(lcsBMatrix,2,'zscore');
%(feltPadMatrix - min(feltPadMatrix, [], 2))./(max(feltPadMatrix, [], 2) - min(feltPadMatrix, [], 2));

[~, idx] = max(normedlcsBMatrix, [], 2); % Returns 2 arrays. Val is the max value of each ROI, idx is the column index of the max value of each ROI.
idxMatrix = horzcat(idx, normedlcsBMatrix); % Appends idx column array to the ROI matrix.
sortedMatrix = sortrows(idxMatrix, 1); % Sorts ROI matrix by the index of their max values in ascending order.
sortedMatrix(: , 1) = []; % Removes indexes from first column of sortedMatrix.
imagesc(sortedMatrix)
set(gca,'clim',[0 3])
xline(42/1.5,'w','LineWidth',1)
xline(47/1.5,'w','LineWidth',1)
xline(122/1.5,'w','LineWidth',1)
xline(127/1.5,'w','LineWidth',1)

xline(63/1.5,'w','LineWidth',1)
xline(68/1.5,'w','LineWidth',1)
xline(100/1.5,'w','LineWidth',1)
xline(106/1.5,'w','LineWidth',1)
ax = gca;
ax.XColor = [0 0 0];
ax.YColor = [0 0 0];
ax.XLabel.Color = [0 0 0];
ax.YLabel.Color = [0 0 0];
set(gca, 'FontSize',18, 'FontName', 'Arial');
ylabel('Cell #')
yticks([1, size(sortedMatrix,1)])
xticks([0 33 67 105])
xticklabels({'0','50','100','157'})

subplot(1,3,2)
lcsAMatrix = cell2mat(lcsA_rmap_mean')';
normedlcsAMatrix   = normalize(lcsAMatrix,2,'zscore');
% normedSandpaperMatrix    = (sandpaperMatrix - min(sandpaperMatrix, [], 2))./(max(sandpaperMatrix, [], 2) - min(sandpaperMatrix, [], 2));

[~, idx] = max(normedlcsAMatrix, [], 2); % Returns 2 arrays. Val is the max value of each ROI, idx is the column index of the max value of each ROI.
idxMatrix = horzcat(idx, normedlcsAMatrix); % Appends idx column array to the ROI matrix.
sortedMatrix = sortrows(idxMatrix, 1); % Sorts ROI matrix by the index of their max values in ascending order.
sortedMatrix(: , 1) = []; % Removes indexes from first column of sortedMatrix.
imagesc(sortedMatrix)

xline(42/1.5,'w','LineWidth',1)
xline(47/1.5,'w','LineWidth',1)
xline(122/1.5,'w','LineWidth',1)
xline(127/1.5,'w','LineWidth',1)

xline(63/1.5,'w','LineWidth',1)
xline(68/1.5,'w','LineWidth',1)
xline(100/1.5,'w','LineWidth',1)
xline(106/1.5,'w','LineWidth',1)
set(gca,'clim',[0 3])
ax = gca;
ax.XColor = [0 0 0];
ax.YColor = [0 0 0];
ax.XLabel.Color = [0 0 0];
ax.YLabel.Color = [0 0 0];
set(gca, 'FontSize',18, 'FontName', 'Arial');
ylabel('Cell #')
yticks([1, size(sortedMatrix,1)])
xticks([1 33 67 105])
xticklabels({'0','50','100','157'})


subplot(1,3,3)
bothMatrix = cell2mat(lcsAB_rmap_mean')';
normedBothMatrix   = normalize(bothMatrix,2,'zscore');

[~, idx] = max(normedBothMatrix, [], 2); % Returns 2 arrays. Val is the max value of each ROI, idx is the column index of the max value of each ROI.
idxMatrix = horzcat(idx, normedBothMatrix); % Appends idx column array to the ROI matrix.
sortedMatrix = sortrows(idxMatrix, 1); % Sorts ROI matrix by the index of their max values in ascending order.
sortedMatrix(: , 1) = []; % Removes indexes from first column of sortedMatrix.
imagesc(sortedMatrix)

xline(42/1.5,'w','LineWidth',1)
xline(47/1.5,'w','LineWidth',1)
xline(122/1.5,'w','LineWidth',1)
xline(127/1.5,'w','LineWidth',1)
xline(63/1.5,'w','LineWidth',1)
xline(68/1.5,'w','LineWidth',1)
xline(100/1.5,'w','LineWidth',1)
xline(106/1.5,'w','LineWidth',1)
set(gca,'clim',[0 3])
ax = gca;
ax.XColor = [0 0 0];
ax.YColor = [0 0 0];
ax.XLabel.Color = [0 0 0];
ax.YLabel.Color = [0 0 0];
set(gca, 'FontSize',18, 'FontName', 'Arial');
ylabel('Cell #')
yticks([1, size(sortedMatrix,1)])
xticks([1 33 67 105])
xticklabels({'0','50','100','157'})


saveas(figSequence,fullfile(savedrct,'sequence.fig'),'fig')
saveas(figSequence,fullfile(savedrct,'sequence.png'),'png')
saveas(figSequence,fullfile(savedrct,'sequence.pdf'),'pdf')


figMean = figure();
subplot(3,1,1)
set(figMean, 'Units', 'centimeters');
set(figMean, 'Position', [0 0 18 8]);
plot(nanmean(normedlcsBMatrix),'k','LineWidth',1.5)
ax = gca;
ax.XColor = [0 0 0];
ax.YColor = [0 0 0];
ax.XLabel.Color = [0 0 0];
ax.YLabel.Color = [0 0 0];
set(gca, 'FontSize',18, 'FontName', 'Arial');
box off
xticks([1 33 67 105])
xticklabels({'0','50','100','157'})
xlim([1 105])
yticks([])

xline(42/1.5,'k','LineWidth',1)
xline(47/1.5,'k','LineWidth',1)
xline(122/1.5,'k','LineWidth',1)
xline(127/1.5,'k','LineWidth',1)

xline(63/1.5,'k','LineWidth',1)
xline(68/1.5,'k','LineWidth',1)
xline(100/1.5,'k','LineWidth',1)
xline(106/1.5,'k','LineWidth',1)

subplot(3,1,2)
plot(nanmean(normedlcsAMatrix),'k','LineWidth',1.5)
ax = gca;
ax.XColor = [0 0 0];
ax.YColor = [0 0 0];
ax.XLabel.Color = [0 0 0];
ax.YLabel.Color = [0 0 0];
set(gca, 'FontSize',18, 'FontName', 'Arial');
box off
xticks([1 33 67 105])
xticklabels({'0','50','100','157'})
xlim([1 105])
yticks([])


xline(42/1.5,'k','LineWidth',1)
xline(47/1.5,'k','LineWidth',1)
xline(122/1.5,'k','LineWidth',1)
xline(127/1.5,'k','LineWidth',1)

xline(63/1.5,'k','LineWidth',1)
xline(68/1.5,'k','LineWidth',1)
xline(100/1.5,'k','LineWidth',1)
xline(106/1.5,'k','LineWidth',1)

subplot(3,1,3)
plot(nanmean(normedBothMatrix),'k','LineWidth',1.5)
ax = gca;
ax.XColor = [0 0 0];
ax.YColor = [0 0 0];
ax.XLabel.Color = [0 0 0];
ax.YLabel.Color = [0 0 0];
set(gca, 'FontSize',18, 'FontName', 'Arial');
box off
xticks([1 33 67 105])
xticklabels({'0','50','100','157'})
xlim([1 105])
yticks([])


xline(42/1.5,'k','LineWidth',1)
xline(47/1.5,'k','LineWidth',1)
xline(122/1.5,'k','LineWidth',1)
xline(127/1.5,'k','LineWidth',1)

xline(63/1.5,'k','LineWidth',1)
xline(68/1.5,'k','LineWidth',1)
xline(100/1.5,'k','LineWidth',1)
xline(106/1.5,'k','LineWidth',1)

saveas(figMean,fullfile(savedrct,'meanCue.fig'),'fig')
saveas(figMean,fullfile(savedrct,'meanCue.png'),'png')
saveas(figMean,fullfile(savedrct,'meanCue.pdf'),'pdf')

