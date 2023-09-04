%% For persistent landmark cells, show the distribution of change in peak firing position when landmark is omitted vs not
[~,peak_ind_off] = max(all_tcs_off(:,1:60)');
[~,peak_ind_on] = max(all_tcs_on(:,1:60)');

peak_change = (peak_ind_off-peak_ind_on)*1.5;
%peak_change = (peak_ind_off-32)*1.5;
figure(4); clf;

histogram(peak_change,'FaceColor',[0.8,.8,.8])
hold on;
line([0,0],[0,60],'Color','r','LineWidth',2)
xlabel('Change in position (cm)')
ylabel('# cells')
title('Peak response position landmark vs no landmark')
set(gca,'FontSize',16)

% subplot(1,2,2);
% histogram((peak_ind_off-32)*1.5)
% hold on
% histogram((peak_ind_on-32)*1.5)

%% Plot sorted peakmaps of all place cells before and after omitting the landmark
all_pcs_on = [];
all_pcs_off = [];
for s = 1:length(mData)
    
    all_pcs_on = [all_pcs_on; sb.generate.averagedTuningFromRasterMaps(mData(s).rmaps.pcs.deconv_motor_on,'smoothing',0)];
    all_pcs_off = [all_pcs_off; sb.generate.averagedTuningFromRasterMaps(mData(s).rmaps.pcs.deconv_motor_off,'smoothing',0)];
    
end

[all_pcs_on,inds] = sb.sort.rasterMapByPeak(all_pcs_on);
all_pcs_off = sb.sort.rasterMapByPeak(all_pcs_off,'sort_order',inds);

% Generate figure
figure(3);
clf;

% Place cells, both landmarks present
subplot(1,2,1);
imagesc(all_pcs_on);
colormap('jet');
line([33,33],[0,2000],'Color','w','LineWidth',line_width)
line([73,73],[0,2000],'Color','w','LineWidth',line_width)
hold on
ylim([0,size(all_pcs_on,1)+250])
plot((size(all_pcs_on,1)+250)-normalize(nanmean(all_pcs_on),'range',[0,200]),'Color','k','LineWidth',line_width)
yticks([0:500:sum(data.n_pcs),sum(data.n_pcs)])

xtix = get(gca, 'XTick');
set(gca, 'XTick',xtix, 'XTickLabel',xtix*1.5)
xticks([0,33,73,105])
xticklabels({'0','50','110','157'})
xlabel('Position (cm)');
title('Stimulus ON');
set(gca,'FontSize',18);

% Place cells, first landmark omitted
subplot(1,2,2);
imagesc(all_pcs_off);
colormap('jet');
line([33,33],[0,2000],'Color','w','LineWidth',line_width,'LineStyle','--')
line([73,73],[0,2000],'Color','w','LineWidth',line_width)
hold on

ylim([0,size(all_pcs_on,1)+250])
plot((size(all_pcs_on,1)+250)-normalize(nanmean(all_pcs_off),'range',[0,200]),'Color','k','LineWidth',line_width)
yticks([]);
xtix = get(gca, 'XTick');
set(gca, 'XTick',xtix, 'XTickLabel',xtix*1.5)
xticks([0,33,73,105])
xticklabels({'0','50','110','157'})
xlabel('Position (cm)');
title('Stimulus OFF');
set(gca,'FontSize',18);

% Set figure size
set(gcf,'renderer', 'painters', 'Position', [1000,1000,450,300])


%% Does the second landmark have a larger response if the first is omitted?
% Get average response of all landmark cells (persistent (tcs) and non-persistent (lcs))
all_tcs_on = [];
all_tcs_off = [];
for s = 1:length(mData)
    all_tcs_on = [all_tcs_on; sb.generate.averagedTuningFromRasterMaps(mData(s).rmaps.tcs.dff_motor_on,'smoothing',0)];
    all_tcs_off = [all_tcs_off; sb.generate.averagedTuningFromRasterMaps(mData(s).rmaps.tcs.dff_motor_off,'smoothing',0)];
end

all_lcs_on = [];
all_lcs_off = [];
for s = 1:length(mData)
    all_lcs_on = [all_lcs_on; sb.generate.averagedTuningFromRasterMaps(mData(s).rmaps.lcs.dff_motor_on,'smoothing',0)];
    all_lcs_off = [all_lcs_off; sb.generate.averagedTuningFromRasterMaps(mData(s).rmaps.lcs.dff_motor_off,'smoothing',0)];
end

% Compare the max response at the second landmark position for ON vs OFF
lcs_peak_on = max(all_lcs_on(:,65:90)');
lcs_peak_off = max(all_lcs_off(:,65:90)');
tcs_peak_on = max(all_tcs_on(:,65:90)');
tcs_peak_off = max(all_tcs_off(:,65:90)');

figure(6); clf; 
colors = cbrewer('div','PRGn',5);
histogram(tcs_peak_on-tcs_peak_off,'BinEdges',[-0.3:0.02:0.3],'FaceColor',colors(1,:),'FaceAlpha',0.8,'Normalization','probability')
hold on
histogram(lcs_peak_on-lcs_peak_off,'BinEdges',[-0.3:0.02:0.3],'FaceColor',colors(5,:),'FaceAlpha',0.7,'Normalization','probability')
legend({'Persistent','Non-persistent'})

yticks([0,0.1,0.2])
ylim([0,0.2])
yticklabels({0, 10, 20})
ylabel('% of cells')
xlabel('? dF/F')
xticks([-0.3,-0.15,0,0.15,0.3])
title('ON vs OFF at second landmark')
set(gca,'FontSize',16);

%% Quantify the amplitude to persistent and non-persistent landmark cells when the landmark is gone, compared to place cells that also have their field close to the landmark
% We need to find all cells that have a firing close to the landmark. We
% should probably use the dff for this:

% Get average response of all landmark cells (persistent (tcs) and non-persistent (lcs))
all_tcs_on = [];
all_tcs_off = [];
for s = 1:length(mData)
    all_tcs_on = [all_tcs_on; sb.generate.averagedTuningFromRasterMaps(mData(s).rmaps.tcs.dff_motor_on,'smoothing',0)];
    all_tcs_off = [all_tcs_off; sb.generate.averagedTuningFromRasterMaps(mData(s).rmaps.tcs.dff_motor_off,'smoothing',0)];
end

all_lcs_on = [];
all_lcs_off = [];
for s = 1:length(mData)
    all_lcs_on = [all_lcs_on; sb.generate.averagedTuningFromRasterMaps(mData(s).rmaps.lcs.dff_motor_on,'smoothing',0)];
    all_lcs_off = [all_lcs_off; sb.generate.averagedTuningFromRasterMaps(mData(s).rmaps.lcs.dff_motor_off,'smoothing',0)];
end

all_pcs_on = [];
all_pcs_off = [];
for s = 1:length(mData)
    all_pcs_on = [all_pcs_on; sb.generate.averagedTuningFromRasterMaps(mData(s).rmaps.pcs.dff_motor_on,'smoothing',0)];
    all_pcs_off = [all_pcs_off; sb.generate.averagedTuningFromRasterMaps(mData(s).rmaps.pcs.dff_motor_off,'smoothing',0)];
end

% Include only place cells with peak activity at landmark
[all_pcs_on_sorted,sort_ind] = sb.sort.rasterMapByPeak(all_pcs_on);
[all_pcs_off_sorted] = sb.sort.rasterMapByPeak(all_pcs_off,'sort_order',sort_ind);

%sb.plot.rasterMaps(cat(3,all_pcs_on_sorted,all_pcs_off_sorted))

% Compare the max response at the second landmark position for ON vs OFF
bins_to_include = 22:62;
lcs_on_mean = nanmax(all_lcs_on(:,bins_to_include)');
lcs_off_mean = nanmax(all_lcs_off(:,bins_to_include)');
tcs_on_mean = nanmax(all_tcs_on(:,bins_to_include)');
tcs_off_mean = nanmax(all_tcs_off(:,bins_to_include)');

pcs_on_mean = nanmax(all_pcs_on(:,bins_to_include)');
pcs_off_mean = nanmax(all_pcs_off(:,bins_to_include)');

figure(7); clf;
grp = [ones(1,length(lcs_on_mean)),ones(1,length(lcs_off_mean))*2,ones(1,length(tcs_on_mean))*3,ones(1,length(tcs_off_mean))*4];
boxplot([lcs_on_mean,lcs_off_mean,tcs_on_mean,tcs_off_mean],grp)
%ylim([-0.001,0.02])

% 
% figure(7); clf;
% b = bar([[nanmean(lcs_on_mean);nanmean(lcs_off_mean)],[nanmean(tcs_on_mean);nanmean(tcs_off_mean)],[nanmean(pcs_on_mean);nanmean(pcs_off_mean)]]');

%% Are place cells affected by the omitted landmark?
all_pcs_on = [];
all_pcs_off = [];
for s = 1:length(mData)
    all_pcs_on = [all_pcs_on; sb.generate.averagedTuningFromRasterMaps(mData(s).rmaps.pcs.deconv_motor_on,'smoothing',0)];
    all_pcs_off = [all_pcs_off; sb.generate.averagedTuningFromRasterMaps(mData(s).rmaps.pcs.deconv_motor_off,'smoothing',0)];
end

% Include only place cells with peak activity at landmark
[all_pcs_on_sorted,sort_ind] = sb.sort.rasterMapByPeak(all_pcs_on);
[all_pcs_off_sorted] = sb.sort.rasterMapByPeak(all_pcs_off,'sort_order',sort_ind);

% Look at the amplitude difference in on vs off 
[~,max_ind] = nanmax(all_pcs_on');
bins_to_include = 22:62;
selected_place_cells_on = all_pcs_on(find(ismember(max_ind,bins_to_include)),:);
selected_place_cells_off = all_pcs_off(find(ismember(max_ind,bins_to_include)),:);

[~,pcs_off_max_ind] = max(selected_place_cells_off');
[~,pcs_on_max_ind] = max(selected_place_cells_on');
pcs_shift = (pcs_off_max_ind-pcs_on_max_ind)*1.5;

figure(1);clf;
subplot(1,2,1);
for i = 1:length(pcs_shift)
    
    if(pcs_shift(i)>0)
        plot([0,pcs_shift(i)],[i,i],'Color','b','LineWidth',3)
    elseif(pcs_shift(i)<0)
        plot([pcs_shift(i),0],[i,i],'Color','r','LineWidth',3)
    else
       plot([-0.5,.5],[i,i],'Color','k','LineWidth',3)
    end
    
    hold on
end

ylim([-length(pcs_shift)*0.05,length(pcs_shift)+(length(pcs_shift)*0.05)]);
yticks([0,length(pcs_shift)]);
yticklabels
xlabel('Shift in position (cm)')
ylabel('Cell #');
title('Place cells')
set(gca,'FontSize',16)

% Landmark cells persistent
all_tcs_on = [];
all_tcs_off = [];
for s = 1:length(mData)
    all_tcs_on = [all_tcs_on; sb.generate.averagedTuningFromRasterMaps(mData(s).rmaps.tcs.deconv_motor_on,'smoothing',0)];
    all_tcs_off = [all_tcs_off; sb.generate.averagedTuningFromRasterMaps(mData(s).rmaps.tcs.deconv_motor_off,'smoothing',0)];
end

% Include only place cells with peak activity at landmark
[all_tcs_on_sorted,sort_ind] = sb.sort.rasterMapByPeak(all_tcs_on);
[all_tcs_off_sorted] = sb.sort.rasterMapByPeak(all_tcs_off,'sort_order',sort_ind);

% Look at the amplitude difference in on vs off 
[~,max_ind] = nanmax(all_tcs_on');
selected_place_cells_on = all_tcs_on(find(ismember(max_ind,bins_to_include)),:);
selected_place_cells_off = all_tcs_off(find(ismember(max_ind,bins_to_include)),:);

[~,tcs_off_max_ind] = max(selected_place_cells_off');
[~,tcs_on_max_ind] = max(selected_place_cells_on');
tcs_shift = (tcs_off_max_ind-tcs_on_max_ind)*1.5;

subplot(1,2,2);
for i = 1:length(tcs_shift)
    
    if(tcs_shift(i)>0)
        plot([0,tcs_shift(i)],[i,i],'Color','b','LineWidth',5)
    elseif(tcs_shift(i)<0)
        plot([tcs_shift(i),0],[i,i],'Color','r','LineWidth',5)
    else
       plot([-0.5,.5],[i,i],'Color','k','LineWidth',5)
    end
    
    hold on
end

ylim([-length(tcs_shift)*0.05,length(tcs_shift)+(length(tcs_shift)*0.05)]);
yticks([0,length(tcs_shift)]);
yticklabels
xlabel('Shift in position (cm)')
ylabel('Cell #');
title('Persistent landmark cells')
set(gca,'FontSize',16)

subplot(1,2,1);
xlim([-110,110])

subplot(1,2,2);
xlim([-110,110])


%% Plot the start and end point of the shift of place and persistent landmark cells when the landmark is omitted

all_pcs_on = [];
all_pcs_off = [];
for s = 1:length(mData)
    all_pcs_on = [all_pcs_on; sb.generate.averagedTuningFromRasterMaps(mData(s).rmaps.pcs.deconv_motor_on,'smoothing',0)];
    all_pcs_off = [all_pcs_off; sb.generate.averagedTuningFromRasterMaps(mData(s).rmaps.pcs.deconv_motor_off,'smoothing',0)];
end

% Include only place cells with peak activity at landmark
[all_pcs_on_sorted,sort_ind] = sb.sort.rasterMapByPeak(all_pcs_on);
[all_pcs_off_sorted] = sb.sort.rasterMapByPeak(all_pcs_off,'sort_order',sort_ind);

% Look at the amplitude difference in on vs off 
[~,max_ind] = nanmax(all_pcs_on');
bins_to_include = 1:105;
selected_place_cells_on = all_pcs_on(find(ismember(max_ind,bins_to_include)),:);
selected_place_cells_off = all_pcs_off(find(ismember(max_ind,bins_to_include)),:);

[~,pcs_off_max_ind] = max(selected_place_cells_off(:,:)');
[~,pcs_on_max_ind] = max(selected_place_cells_on(:,:)');
pcs_shift = (pcs_off_max_ind-pcs_on_max_ind)*1.5;

% Sort them
[pcs_shift,sort_ind] = sort(pcs_shift);
pcs_on_max_ind = pcs_on_max_ind(sort_ind);
pcs_off_max_ind = pcs_off_max_ind(sort_ind);

figure(1);clf;
subplot(1,2,1);
for i = 1:length(pcs_shift)
    if(pcs_shift(i)>0)
        plot([pcs_on_max_ind(i),pcs_off_max_ind(i)],[i,i],'Color','b','LineWidth',3)
    elseif(pcs_shift(i)==0)
        plot([pcs_on_max_ind(i)-0.5,pcs_off_max_ind(i)+0.5],[i,i],'Color','k','LineWidth',3)
    else
        plot([pcs_on_max_ind(i),pcs_off_max_ind(i)],[i,i],'Color','r','LineWidth',3)
    end
    hold on
end

ylim([-length(pcs_shift)*0.05,length(pcs_shift)+(length(pcs_shift)*0.05)]);
yticks([0,length(pcs_shift)]);
yticklabels
xlabel('Shift in position (cm)')
ylabel('Cell #');
title('Place cells')
xlim([0,105])
set(gca,'FontSize',16)


all_tcs_on = [];
all_tcs_off = [];
for s = 1:length(mData)
    all_tcs_on = [all_tcs_on; sb.generate.averagedTuningFromRasterMaps(mData(s).rmaps.tcs.deconv_motor_on,'smoothing',0)];
    all_tcs_off = [all_tcs_off; sb.generate.averagedTuningFromRasterMaps(mData(s).rmaps.tcs.deconv_motor_off,'smoothing',0)];
end

% Include only place cells with peak activity at landmark
[all_tcs_on_sorted,sort_ind] = sb.sort.rasterMapByPeak(all_tcs_on);
[all_tcs_off_sorted] = sb.sort.rasterMapByPeak(all_tcs_off,'sort_order',sort_ind);

% Look at the amplitude difference in on vs off 
[~,max_ind] = nanmax(all_tcs_on');
bins_to_include = 1:105;
selected_place_cells_on = all_tcs_on(find(ismember(max_ind,bins_to_include)),:);
selected_place_cells_off = all_tcs_off(find(ismember(max_ind,bins_to_include)),:);

[~,tcs_off_max_ind] = max(selected_place_cells_off(:,bins_to_include)');
[~,tcs_on_max_ind] = max(selected_place_cells_on(:,bins_to_include)');
tcs_shift = (tcs_off_max_ind-tcs_on_max_ind);

% Sort them
[tcs_shift,sort_ind] = sort(tcs_shift);
tcs_on_max_ind = tcs_on_max_ind(sort_ind);
tcs_off_max_ind = tcs_off_max_ind(sort_ind);

subplot(1,2,2);
for i = 1:length(tcs_shift)
    if(tcs_shift(i)>0)
        plot([tcs_on_max_ind(i),tcs_off_max_ind(i)],[i,i],'Color','b','LineWidth',3)
    elseif(tcs_shift(i)==0)
        plot([tcs_on_max_ind(i)-0.5,tcs_off_max_ind(i)+0.5],[i,i],'Color','b','LineWidth',3)
    else
        plot([tcs_on_max_ind(i),tcs_off_max_ind(i)],[i,i],'Color','r','LineWidth',3)
    end
    hold on
end

ylim([-length(tcs_shift)*0.05,length(tcs_shift)+(length(tcs_shift)*0.05)]);
yticks([0,length(tcs_shift)]);
yticklabels
xlabel('Shift in position (cm)')
ylabel('Cell #');
title('Landmark cells')
xlim([0,105])
set(gca,'FontSize',16)



%% Is there a drop in % place cells in the OFF laps vs the ON laps?
% In order to test this, I use odd vs even laps in ON laps as a control.

% Compute the percentage of place cells 
pcs_on_odd = [];
pcs_on_even = [];
pcs_off = [];
for s = 1:length(mData)
    
    pcs_on_odd(s) = length(sb.classify.placeTunedRoisSimple(mData(s).rmaps.deconv_motor_on(1:2:end,:,:)))/size(mData(s).rmaps.dff,3);
    pcs_on_even(s) = length(sb.classify.placeTunedRoisSimple(mData(s).rmaps.deconv_motor_on(2:2:end,:,:)))/size(mData(s).rmaps.dff,3);
    pcs_off(s) = length(sb.classify.placeTunedRoisSimple(mData(s).rmaps.deconv_motor_off))/size(mData(s).rmaps.dff,3);
    
end

figure(1); clf;
plot(pcs_on_odd)
hold on;
plot(pcs_on_even);
plot(pcs_off);
ylim([0,0.60])
title('% place cells in each session');
legend({'Odd','Even','Off'})
xlabel('Session #');
xlim([0,11])
set(gca,'FontSize',16);


%% Do place cells loose their place cell tuning when the landmark is omitted?
pcs_diff_odd_vs_even = [];
pcs_diff_odd_vs_off = [];
for s = 1:length(mData)
    
    pcs_odd = sb.classify.placeTunedRoisSimple(mData(s).rmaps.deconv_motor_on(1:2:end,:,:));
    pcs_even = sb.classify.placeTunedRoisSimple(mData(s).rmaps.deconv_motor_on(2:2:end,:,:));
    pcs_off = sb.classify.placeTunedRoisSimple(mData(s).rmaps.deconv_motor_off(:,:,:));
    
    pcs_diff_odd_vs_even(s) = sum(ismember(pcs_odd,pcs_even))/length(pcs_odd);
    pcs_diff_odd_vs_off(s) = sum(ismember(pcs_odd,pcs_off))/length(pcs_odd);
    
end

figure(1);clf;
bar([nanmean(pcs_diff_odd_vs_even),nanmean(pcs_diff_odd_vs_off)])
hold on;

% Add errorbars
sem = std(pcs_diff_odd_vs_even)/sqrt(length(pcs_diff_odd_vs_even));
er = errorbar(1,nanmean(pcs_diff_odd_vs_even),-sem,sem);
er.Color = [0,0,0];
er.LineStyle = 'none';
er.LineWidth = 2;

sem = std(pcs_diff_odd_vs_off)/sqrt(length(pcs_diff_odd_vs_off));
er = errorbar(2,nanmean(pcs_diff_odd_vs_off),-sem,sem);
er.Color = [0,0,0];
er.LineStyle = 'none';
er.LineWidth = 2;
ylim([0,1])

xticklabels({'Control', 'On vs Off'})

set(gca,'FontSize',16)


%% Place cells - Compute the change in amplitude 

all_pcs_on = [];
all_pcs_off = [];
for s = 1:length(mData)
    all_pcs_on = [all_pcs_on; sb.generate.averagedTuningFromRasterMaps(mData(s).rmaps.pcs.deconv_motor_on,'smoothing',0)];
    all_pcs_off = [all_pcs_off; sb.generate.averagedTuningFromRasterMaps(mData(s).rmaps.pcs.deconv_motor_off,'smoothing',0)];
end



figure(1); clf;
plot(nanmean(all_pcs_on))
hold on
plot(nanmean(all_pcs_off))






