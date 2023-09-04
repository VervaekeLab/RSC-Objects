
% How many rois are included in total?
n_acs = zeros(1,3);

for s = 1:length(mData)
    
    n_acs(1) = n_acs(1) + length(mData(s).rmaps.level1.acs);
    n_acs(2) = n_acs(2) + length(mData(s).rmaps.level2.acs);
    n_acs(3) = n_acs(3) + length(mData(s).rmaps.level3.acs);
    
end


%% ------ LANDMARK CELLS
%% Show the percentage of cells that respond in each category
p_lcs_1 = [];
p_lcs_2 = [];
p_lcs_3 = [];
for s = 1:length(mData)
    p_lcs_1(s) = (length(mData(s).rmaps.level1.lcs) / size(mData(s).rmaps.level1.dff,3));
    p_lcs_2(s) = (length(mData(s).rmaps.level2.lcs) / size(mData(s).rmaps.level2.dff,3));
    p_lcs_3(s) = (length(mData(s).rmaps.level3.lcs) / size(mData(s).rmaps.level3.dff,3));
end

p_lcs_1 = p_lcs_1(1:end);
p_lcs_2 = p_lcs_2(1:end);
p_lcs_3 = p_lcs_3(1:end);

% Plot bar plot with error
x = 1:3;
data = [nanmean(p_lcs_1),nanmean(p_lcs_2),nanmean(p_lcs_3)];
errhigh = [nanstd(p_lcs_1)/sqrt(length(p_lcs_1)), nanstd(p_lcs_2)/sqrt(length(p_lcs_2)), nanstd(p_lcs_3)/sqrt(length(p_lcs_3))];
errlow = errhigh;

figure(3);
clf;
subplot(1,3,1);
landmark_cells_color = [239,35,60]/255;
bar(x,data,'FaceColor',landmark_cells_color)

hold on
er = errorbar(x,data,errlow,errhigh);
er.Color = [0 0 0];
er.LineStyle = 'none';

yticks([0.05, 0.1,0.15,0.2])
yticklabels({5,10,15,20})
xticks([1,2,3])
title('Tactile cells')
ylim([0,0.24])
xtickangle(45)
xticklabels({'No light','Dim light','Bright light'})
ylabel('% of cells')
set(gca,'FontSize',16)


%% ------ PLACE CELLS
%% Show the percentage of cells that respond in each category
p_pcs_1 = [];
p_pcs_2 = [];
p_pcs_3 = [];
for s = 1:length(mData)
    p_pcs_1(s) = (length(mData(s).rmaps.level1.pcs) / size(mData(s).rmaps.level1.dff,3));
    p_pcs_2(s) = (length(mData(s).rmaps.level2.pcs) / size(mData(s).rmaps.level2.dff,3));
    p_pcs_3(s) = (length(mData(s).rmaps.level3.pcs) / size(mData(s).rmaps.level3.dff,3));   
end

p_pcs_1 = p_pcs_1(1:end);
p_pcs_2 = p_pcs_2(1:end);
p_pcs_3 = p_pcs_3(1:end);

% Plot bar plot with error
x = 1:3;
data = [nanmean(p_pcs_1),nanmean(p_pcs_2),nanmean(p_pcs_3)];
errhigh = [nanstd(p_pcs_1)/sqrt(length(p_pcs_1)), nanstd(p_pcs_2)/sqrt(length(p_pcs_2)), nanstd(p_pcs_3)/sqrt(length(p_pcs_3))];
errlow = errhigh;

place_cells_color = [51,51,51]/255;
subplot(1,3,2);
bar(x,data,'FaceColor',place_cells_color)

hold on
er = errorbar(x,data,errlow,errhigh);
er.Color = [0 0 0];
er.LineStyle = 'none';

yticks([0.1, 0.2,0.3])
yticklabels({10,20,30})
xticks([1,2,3])
title('Position-tuned cells')
ylim([0,0.33])
xtickangle(45)
xticklabels({'No light','Dim light','Bright light'})
ylabel('% of cells')
set(gca,'FontSize',16)

% Set correct size of the figure
set(gcf,'renderer', 'painters', 'Position', [1,100,800,800])



%% Plot the change in the number of place cells and landmark cells compared to the baseline of no light.
n_pcs = [];
n_lcs = [];
for s = 1:length(mData)
    levels = {'level1','level2','level3'};
   for l = 1:3
        n_pcs(s,l) = length(mData(s).rmaps.(levels{l}).pcs) /size(mData(s).rmaps.(levels{l}).dff,3);
        n_lcs(s,l) = length(mData(s).rmaps.(levels{l}).lcs) /size(mData(s).rmaps.(levels{l}).dff,3);
   end
end

% Convert to percent
n_pcs = n_pcs*100;
n_lcs = n_lcs*100;

landmark_cells_color = [239,35,60]/255;
place_cells_color = [51,51,51]/255;

subplot(1,3,3);
n_lcs_mean = nanmean(n_lcs);
change_in_lcs = (n_lcs_mean./n_lcs_mean(1))*100;
n_pcs_mean = nanmean(n_pcs);
change_in_pcs = (n_pcs_mean./n_pcs_mean(1))*100;
y = [change_in_lcs;change_in_pcs]';
b = bar(y);
b(1).FaceColor = landmark_cells_color;
b(2).FaceColor = place_cells_color;
yticks([0:100:600])
ylim([0,350])
xticklabels({'No light','Dim light','Bright light'})
xtickangle(45)
ylabel('% of cells')
title('Change in number of cells')
l = legend({'Tactile cells','Position-tuned cells'},'location','northwest');
set(gca,'FontSize',16)

% Add error bars
err_lcs = (nanstd(n_lcs./n_lcs(:,1))/sqrt(5)) * 100;
err_pcs = (nanstd(n_pcs./n_pcs(:,1))/sqrt(5)) * 100;

hold on
errorbar([1,2,3]-0.124,change_in_lcs,err_lcs, '.','Color','k')
errorbar([1,2,3]+0.16,change_in_pcs,err_pcs, '.','Color','k')

l.String = {'Tactile cells','Position-tuned cells'};

% Set correct size of the figure
set(gcf,'renderer', 'painters', 'Position', [100,100,1000,300])

