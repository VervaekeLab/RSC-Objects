%% Create a figure to show the averaged firing of all "landmark cells" and "trace cells" on laps with and without the motor landmark present.
close all;
line_width = 1;
max_z_value = 3;

%% Plot sorted peakmaps of all trace cells before and after omitting the landmark
all_tcs_on = [];
all_tcs_off = [];
for s = 1:length(mData)
    all_tcs_on = [all_tcs_on; sb.generate.averagedTuningFromRasterMaps(mData(s).rmaps.tcs.deconv_motor_on,'smoothing',0)];
    all_tcs_off = [all_tcs_off; sb.generate.averagedTuningFromRasterMaps(mData(s).rmaps.tcs.deconv_motor_off,'smoothing',0)];
end

[all_tcs_on,inds] = sb.sort.rasterMapByPeak(all_tcs_on,'normalize_type','zscore');
all_tcs_off = sb.sort.rasterMapByPeak(all_tcs_off,'sort_order',inds, 'normalize_type','zscore');

% Z score normalize the average firing maps
% all_tcs_on = zScoreNormalize(all_tcs_on,'row');
% all_tcs_off = zScoreNormalize(all_tcs_off,'row');

n_cells = size(all_tcs_on,1)+35; % +35 is added to create a white space at the bottom of the plots where the averaged is plotted

% Generate figure
figure(1);
clf;

% Persistent landmark cells, both landmarks is present
subplot(2,2,3);
imagesc(all_tcs_on,[0,max_z_value]);
line([33,33],[0,1000],'Color','w','LineWidth',line_width)
line([73,73],[0,1000],'Color','w','LineWidth',line_width)
hold on
ylim([0,n_cells])
plot(n_cells-normalize(nanmean(all_tcs_on),'range',[0,32]),'Color','k','LineWidth',line_width)
yticks([1,sum(data.n_tcs)])
ylabel('Cell #');

xtix = get(gca, 'XTick');
set(gca, 'XTick',xtix, 'XTickLabel',xtix*1.5)
xticks([0,33,73,105])
xticklabels({'0','50','110','157'})
xlabel('Position (cm)');
%title('Stimulus ON');
set(gca,'FontSize',18);
ylim([0.5,n_cells])

% Persistent landmark cells, first landmarks is omitted
subplot(2,2,4);
imagesc(all_tcs_off,[0,max_z_value]);
line([33,33],[0,1000],'Color','w','LineWidth',line_width,'LineStyle','--')
line([73,73],[0,1000],'Color','w','LineWidth',line_width)
hold on

ylim([0,n_cells])
plot(n_cells-normalize(nanmean(all_tcs_off),'range',[0,32]),'Color','k','LineWidth',line_width)
yticks([]);
xtix = get(gca, 'XTick');
set(gca, 'XTick',xtix, 'XTickLabel',xtix*1.5)
xticks([0,33,73,105])
xticklabels({'0','50','110','157'})
xlabel('Position (cm)');
%title('Stimulus OFF');
ylim([0.5,n_cells])
set(gca,'FontSize',18);

% Plot sorted peakmaps of all landmark cells before and after omitting the landmark
all_lcs_on = [];
all_lcs_off = [];
no_pcs_all=0;
for s = 1:length(mData)
    
    all_lcs_on = [all_lcs_on; sb.generate.averagedTuningFromRasterMaps(mData(s).rmaps.lcs.deconv_motor_on,'smoothing',0)];
    all_lcs_off = [all_lcs_off; sb.generate.averagedTuningFromRasterMaps(mData(s).rmaps.lcs.deconv_motor_off,'smoothing',0)];
    
    pcs = sb.classify.placeCellsShuffle(mData(s).rmaps.lcs.deconv_motor_off);
    no_pcs_all = no_pcs_all+length(pcs);
end

[all_lcs_on,inds] = sb.sort.rasterMapByPeak(all_lcs_on,'normalize',false);
all_lcs_off = sb.sort.rasterMapByPeak(all_lcs_off,'sort_order',inds,'normalize',false);

all_lcs_on = normalize(all_lcs_on,2);%zScoreNormalize(all_lcs_on,'row');
all_lcs_off = normalize(all_lcs_off,2);%zScoreNormalize(all_lcs_off,'row');

% Non-persistent landmark cells, both landmarks is present
subplot(2,2,1);
imagesc(all_lcs_on,[0,max_z_value]);
line([33,33],[0,1000],'Color','w','LineWidth',line_width)
line([73,73],[0,1000],'Color','w','LineWidth',line_width)
hold on
ylim([0,size(all_lcs_on,1)+25])
plot((size(all_lcs_on,1)+25)-normalize(nanmean(all_lcs_on),'range',[0,20]),'Color','k','LineWidth',line_width)
yticks([1,sum(data.n_lcs)])
ylabel('Cell #');

xtix = get(gca, 'XTick');
set(gca, 'XTick',xtix, 'XTickLabel',xtix*1.5)
xticks([0,33,73,105])
xticklabels({'0','50','110','157'})
%xlabel('Position (cm)');
title('Stimulus ON');
set(gca,'FontSize',18);

% Non-persistent landmark cells, first landmark is omitted
subplot(2,2,2);
imagesc(all_lcs_off,[0,max_z_value]);
line([33,33],[0,1000],'Color','w','LineWidth',line_width,'LineStyle','--')
line([73,73],[0,1000],'Color','w','LineWidth',line_width)
hold on
ylim([0,size(all_lcs_on,1)+25])
plot((size(all_lcs_on,1)+25)-normalize(nanmean(all_lcs_off),'range',[0,20]),'Color','k','LineWidth',line_width)
yticks([]);
xtix = get(gca, 'XTick');
set(gca, 'XTick',xtix, 'XTickLabel',xtix*1.5)
xticks([0,33,73,105])
xticklabels({'0','50','110','157'})
%xlabel('Position (cm)');
title('Stimulus OFF');
set(gca,'FontSize',18);

% Set figure size
set(gcf,'renderer', 'painters', 'Position', [-1200,1000,550,800])



