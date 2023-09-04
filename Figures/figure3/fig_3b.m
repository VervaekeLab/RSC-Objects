% Example of landmark cells that are exclusively active in dark or
% exclusively in light

figure(2); clf;
signal_type = 'deconv';
line_width = 0.5;
z_score_level = 4;

%  --- Example cell 1
subplot(1,3,1);
s = 3;
c = 28; % this is cell 163
rmap1 = mData(s).rmaps.level1.(signal_type)(:,:,mData(s).rmaps.level1.lcs(c));
rmap2 = mData(s).rmaps.level2.(signal_type)(:,:,mData(s).rmaps.level1.lcs(c));
rmap3 = mData(s).rmaps.level3.(signal_type)(:,:,mData(s).rmaps.level1.lcs(c));

rmap = cat(1,rmap1,rmap2,rmap3);
rmap = smoothdata(rmap,2,'gaussian',5);

% Zscore normalize
% rmap = zScoreNormalize(rmap,'all');
rmap = normalize(rmap);
% Plot
imagesc(rmap,[0,z_score_level])
hold on
colormap('jet')

% Draw lines to indicate change in condition
line([0,105],[51,51],'Color','k')
line([0,105],[101,101],'Color','k')

% Draw lines to indcate landmark
line([32,32],[0,1000],'Color','r');
line([72,72],[0,1000],'Color','r');

% Draw the averages on top
avg_signal1 = normalize(nanmean(rmap1),'range',[0,10]);
avg_signal2 = normalize(nanmean(rmap2),'range',[0,10]);
avg_signal3 = normalize(nanmean(rmap3),'range',[0,10]);
max_max = max([max(nanmean(rmap1)),max(nanmean(rmap2)),max(nanmean(rmap3))]);

% Normalize response to each other
avg_signal1 = avg_signal1 * (max(nanmean(rmap1))/max_max);
avg_signal2 = avg_signal2 * (max(nanmean(rmap2))/max_max);
avg_signal3 = avg_signal3 * (max(nanmean(rmap3))/max_max);

plot(51-avg_signal1,'Color','k','LineWidth',line_width)
plot(101-avg_signal2,'Color','k','LineWidth',line_width)
plot(size(rmap,1)-avg_signal3,'Color','k','LineWidth',line_width)

xticks([1,105])
xticklabels({1,157})
xlabel('Position (cm)')
ylabel('Lap #')
yticks([1,size(rmap,1)])

title('Cell #1')
set(gca,'FontSize',16)

%  --- Example cell 2
subplot(1,3,2);
s = 3;
c = 34; % this is cell 188

rmap1 = mData(s).rmaps.level1.(signal_type)(:,:,188);
rmap2 = mData(s).rmaps.level2.(signal_type)(:,:,188);
rmap3 = mData(s).rmaps.level3.(signal_type)(:,:,188);

rmap1 = mData(s).rmaps.level1.(signal_type)(:,:,mData(s).rmaps.level1.lcs(c));
rmap2 = mData(s).rmaps.level2.(signal_type)(:,:,mData(s).rmaps.level1.lcs(c));
rmap3 = mData(s).rmaps.level3.(signal_type)(:,:,mData(s).rmaps.level1.lcs(c));

rmap = cat(1,rmap1,rmap2,rmap3);
rmap = smoothdata(rmap,2,'gaussian',5);

% Zscore normalize
% rmap = zScoreNormalize(rmap,'all');
rmap = normalize(rmap);
% Plot
imagesc(rmap,[0,z_score_level])
hold on
colormap('jet')

% Draw lines to indicate change in condition
line([0,105],[51,51],'Color','k')
line([0,105],[101,101],'Color','k')

% Draw lines to indcate landmark
line([32,32],[0,1000],'Color','r');
line([72,72],[0,1000],'Color','r');

% Draw the averages on top
avg_signal1 = normalize(nanmean(rmap1),'range',[0,10]);
avg_signal2 = normalize(nanmean(rmap2),'range',[0,10]);
avg_signal3 = normalize(nanmean(rmap3),'range',[0,10]);
max_max = max([max(nanmean(rmap1)),max(nanmean(rmap2)),max(nanmean(rmap3))]);

% Normalize response to each other
avg_signal1 = avg_signal1 * (max(nanmean(rmap1))/max_max);
avg_signal2 = avg_signal2 * (max(nanmean(rmap2))/max_max);
avg_signal3 = avg_signal3 * (max(nanmean(rmap3))/max_max);

plot(51-avg_signal1,'Color','k','LineWidth',line_width)
plot(101-avg_signal2,'Color','k','LineWidth',line_width)
plot(size(rmap,1)-avg_signal3,'Color','k','LineWidth',line_width)


xticks([1,105])
xticklabels({1,157})
xlabel('Position (cm)')
ylabel('Lap #')
yticks([1,size(rmap,1)])
title('Cell #2')
set(gca,'FontSize',16)

%  --- Example cell 3
subplot(1,3,3);
s = 3;
c = 21; % this is cell 116
rmap1 = mData(s).rmaps.level1.(signal_type)(:,:,mData(s).rmaps.level1.lcs(c));
rmap2 = mData(s).rmaps.level2.(signal_type)(:,:,mData(s).rmaps.level1.lcs(c));
rmap3 = mData(s).rmaps.level3.(signal_type)(:,:,mData(s).rmaps.level1.lcs(c));

rmap = cat(1,rmap1,rmap2,rmap3);
rmap = smoothdata(rmap,2,'gaussian',5);

% Zscore normalize
% rmap = zScoreNormalize(rmap,'all');
rmap = normalize(rmap);

% Plot
imagesc(rmap,[0,z_score_level])
hold on
colormap(flipud(colormap('gray')))

% Draw lines to indicate change in condition
line([0,105],[51,51],'Color','k')
line([0,105],[101,101],'Color','k')

% Draw lines to indcate landmark
line([32,32],[0,1000],'Color','r');
line([72,72],[0,1000],'Color','r');

% Draw the averages on top
avg_signal1 = normalize(nanmean(rmap1),'range',[0,10]);
avg_signal2 = normalize(nanmean(rmap2),'range',[0,10]);
avg_signal3 = normalize(nanmean(rmap3),'range',[0,10]);
max_max = max([max(nanmean(rmap1)),max(nanmean(rmap2)),max(nanmean(rmap3))]);

% Normalize response to each other
avg_signal1 = avg_signal1 * (max(nanmean(rmap1))/max_max);
avg_signal2 = avg_signal2 * (max(nanmean(rmap2))/max_max);
avg_signal3 = avg_signal3 * (max(nanmean(rmap3))/max_max);

plot(51-avg_signal1,'Color','k','LineWidth',line_width)
plot(101-avg_signal2,'Color','k','LineWidth',line_width)
plot(size(rmap,1)-avg_signal3,'Color','k','LineWidth',line_width)

xticks([1,105])
xticklabels({1,157})
xlabel('Position (cm)')
ylabel('Lap #')
yticks([1,size(rmap,1)])
title('Cell #3')
set(gca,'FontSize',16)


% Set correct size of the figure
set(gcf,'renderer', 'painters', 'Position', [2500,100,800,400])


%% 



