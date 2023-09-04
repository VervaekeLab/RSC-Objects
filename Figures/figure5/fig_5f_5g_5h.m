%% This creates a plot showing the distribution of place fields around the landmarks.
% We can see that around the landmark, there is a bias towards place fields
% having their peak after the landmark rather than before. Also, there is a
% significantly smaller field width after the landmark compared to before.
% In the following code, "baseline" refers to cells before cue and
% "response" to cells after the cue.

% Place cells are already detected, and I use the first 50 laps where the 
% landmark is still present as data.
place_cells_rmaps_on = [];
place_cells_rmaps_off = [];
n_pcs = 0;
session_ind = [];
for s = 1:length(mData)
   place_cells_rmaps_on = cat(3,place_cells_rmaps_on,mData(s).rmaps.pcs.deconv_motor_on(1:50,:,:));
   place_cells_rmaps_off = cat(3,place_cells_rmaps_off,mData(s).rmaps.pcs.deconv_motor_off(1:27,:,:));
   n_pcs = n_pcs + length(mData(s).rmaps.pcs.ind);
   session_ind = [session_ind,ones(1,length(mData(s).rmaps.pcs.ind))*s];
end

% Get stats about the place cell rastermaps
[~,stats_rmaps_on] = sb.classify.placeTunedRoisSimple(place_cells_rmaps_on);
[~,stats_rmaps_off] = sb.classify.placeTunedRoisSimple(place_cells_rmaps_off);

% Remove all multi peak cells from stats struct.
single_peak_cells = find([stats_rmaps_on.n_peaks] == 1);
stats_rmaps_on = stats_rmaps_on(single_peak_cells);
stats_rmaps_off = stats_rmaps_off(single_peak_cells);
place_cells_rmaps_on = place_cells_rmaps_on(:,:,single_peak_cells);
place_cells_rmaps_off = place_cells_rmaps_off(:,:,single_peak_cells); 

session_ind = session_ind(single_peak_cells);
% Find the peak position of all cells
max_ind = [stats_rmaps_on(:).field_peak_ind];

% Define the window to search
% 20 cm before vs 20 cm after
inds_for_first_baseline = 19:32;
inds_for_first_resp = 33:46;
inds_for_sec_baseline = 59:72;
inds_for_sec_resp = 73:86;

% Plot landmark clustered place cells combined
% Here I combined the place cells at both fields and align the data to the
% lanmdark onset a position 0. 
% Average the activity across all laps for each cells
average_activity_on = sb.generate.averagedTuningFromRasterMaps(place_cells_rmaps_on);
clustered_cells_1 = normalize(average_activity_on(find(ismember(max_ind,19:46)),:)','range',[0,1])';
clustered_cells_2 = normalize(average_activity_on(find(ismember(max_ind,59:86)),:)','range',[0,1])';
clustered_cells = cat(1,clustered_cells_1(:,10:55),clustered_cells_2(:,50:95));
combined_map = sb.sort.rasterMapByPeak(clustered_cells,'normalize',false);

figure(6481);
clf;
subplot(2,2,[1,3])
imagesc(combined_map);
colormap('parula')
line([23.5,23.5],[0,10000],'LineWidth',3,'Color','w')
xticks([10.5,23.5,36.5])
xticklabels({'-20','0','20'})
ylabel('Cell #');
xlabel('Distance from landmark onset (cm)');
title('Place cells close to landmarks');
set(gca,'FontSize',20)

% ------ Find percentage of cells before and after landmark
% Find the number of place cells with peak firing within window
first_baseline = sum(ismember(max_ind,inds_for_first_baseline));
first_resp = sum(ismember(max_ind,inds_for_first_resp));
first_comb = first_baseline+first_resp;
sec_baseline = sum(ismember(max_ind,inds_for_sec_baseline));
sec_resp = sum(ismember(max_ind,inds_for_sec_resp));
sec_comb = sec_baseline+sec_resp;

% Combined first and second landmark
all_baseline = sum(ismember(max_ind,[inds_for_first_baseline,inds_for_sec_baseline]));
all_resp = sum(ismember(max_ind,[inds_for_first_resp,inds_for_sec_resp]));
all_comb = all_baseline + all_resp;
all_baseline = all_baseline / all_comb;
all_resp = all_resp / all_comb;

% Divide the numbers by the sum to get the percentage of cells
first_baseline = first_baseline / first_comb;
first_resp = first_resp / first_comb;

sec_baseline = sec_baseline / sec_comb;
sec_resp = sec_resp / sec_comb;

% Plot figure
subplot(2,2,2)
pie([all_baseline,all_resp])

% ---- Do place cells close to the landmarks have narrower tuning?
first_baseline_inds = find(ismember([stats_rmaps_on(:).field_peak_ind],inds_for_first_baseline));
first_resp_inds = find(ismember([stats_rmaps_on(:).field_peak_ind],inds_for_first_resp));

first_baseline_width_on = [stats_rmaps_on(first_baseline_inds).field_width];
first_resp_width_on = [stats_rmaps_on(first_resp_inds).field_width];

first_baseline_average_width_on = mean(first_baseline_width_on) * 1.5;
first_resp_average_width_on = mean(first_resp_width_on) * 1.5;

% Second landmark
sec_baseline_inds = find(ismember([stats_rmaps_on(:).field_peak_ind],inds_for_sec_baseline));
sec_resp_inds = find(ismember([stats_rmaps_on(:).field_peak_ind],inds_for_sec_resp));

sec_baseline_width = [stats_rmaps_on(sec_baseline_inds).field_width];
sec_resp_width = [stats_rmaps_on(sec_resp_inds).field_width];

sec_baseline_average_width = mean(sec_baseline_width) * 1.5;
sec_resp_average_width = mean(sec_resp_width) * 1.5;

% Lets combine the two landmark positions
baseline_width_on = [first_baseline_width_on, sec_baseline_width];
resp_width_on = [first_resp_width_on,sec_resp_width];

% Plot averages and error
subplot(2,2,4)

baseline_width_on_mean = nanmean(baseline_width_on);
baseline_width_on_sem = nanstd(baseline_width_on)/sqrt(length(baseline_width_on));
resp_width_on_mean = nanmean(resp_width_on);
resp_width_on_sem = nanstd(resp_width_on)/sqrt(length(resp_width_on));

bar([baseline_width_on_mean,resp_width_on_mean])
hold on;
errorbar([1,2],[baseline_width_on_mean,resp_width_on_mean],[baseline_width_on_sem,resp_width_on_sem],'LineStyle','none','Color','k')
ylim([0,30]);
yticks([0,30])
ylabel('Place field width (cm)');
xticklabels({'Before','After'})
title('Place field widths')
set(gca,'FontSize',20);

%  Make the same plot but for laps where the landmark is omitted
% Lets combine the two landmark positions
baseline_width_on = [first_baseline_width_on, sec_baseline_width];
resp_width_on = [first_resp_width_on,sec_resp_width];

set(gcf,'renderer', 'painters', 'Position', [1200 200 1100 1000])

%% Run statistical analysis on field width using Mann-Whitney U test (ranksum)
p = ranksum(resp_width_on,baseline_width_on)



