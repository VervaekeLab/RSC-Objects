% Analysis of place cells during landmark manipulation. Do PCS change their
% field position, reliability, or field width during the laps of landmark
% manipulation?

% What percentage of place cells loose their tuning when the landmark is
% gone? (Can I use this as a proxy to find landmark response cells?)
signal_type_on = "deconv_motor_on";
signal_type_off = "deconv_motor_off";

for s = 1:length(mData)
    
    if s == 1
        pcs_avg_dff_on_odd = sb.generate.averagedTuningFromRasterMaps(mData(s).rmaps.pcs.(signal_type_on)(1:2:end,:,:));
        pcs_avg_dff_on_even = sb.generate.averagedTuningFromRasterMaps(mData(s).rmaps.pcs.(signal_type_on)(2:2:end,:,:));
        pcs_avg_dff_off = sb.generate.averagedTuningFromRasterMaps(mData(s).rmaps.pcs.(signal_type_off));
        pcs_avg_dff_on = sb.generate.averagedTuningFromRasterMaps(mData(s).rmaps.pcs.(signal_type_on));
        
    else
        pcs_avg_dff_on_odd = cat(1,pcs_avg_dff_on_odd,sb.generate.averagedTuningFromRasterMaps(mData(s).rmaps.pcs.(signal_type_on)(1:2:end,:,:)));
        pcs_avg_dff_on_even = cat(1,pcs_avg_dff_on_even,sb.generate.averagedTuningFromRasterMaps(mData(s).rmaps.pcs.(signal_type_on)(2:2:end,:,:)));
        pcs_avg_dff_off = cat(1,pcs_avg_dff_off,sb.generate.averagedTuningFromRasterMaps(mData(s).rmaps.pcs.(signal_type_off)));
        pcs_avg_dff_on = cat(1,pcs_avg_dff_on,sb.generate.averagedTuningFromRasterMaps(mData(s).rmaps.pcs.(signal_type_on)));

    end
  
end



% For each cell, find the amplitude loss around its peak position
[pcs_sorted_on_odd,sort_ind] = sb.sort.rasterMapByPeak(pcs_avg_dff_on_odd,'normalize',true);
[pcs_sorted_on_even] = sb.sort.rasterMapByPeak(pcs_avg_dff_on_even,'sort_order',sort_ind,'normalize',true);
pcs_sorted_off = sb.sort.rasterMapByPeak(pcs_avg_dff_off,'sort_order',sort_ind,'normalize',true);
%sb.plot.rasterMaps(cat(3,pcs_sorted_on_odd,pcs_sorted_on_even,pcs_sorted_off))


%% Plot correlation
figure(1);clf
subplot(1,2,1);
imagesc(corr(pcs_sorted_on_odd,pcs_sorted_on_even));
hold on
line([32,32],[0,2000],'Color','w')
line([72,72],[0,2000],'Color','w')
line([0,2000],[32,32],'Color','w')
line([0,2000],[72,72],'Color','w')
ylabel('Odd')

xlabel('Even');
title('Odd vs Even');
yticks([1,32,72,105]);
yticklabels({'0','50','110','157'});
xticks([1,32,72,105]);
xticklabels({'0','50','110','157'});
set(gca,'FontSize',16)

subplot(1,2,2);
imagesc(corr(pcs_sorted_on_odd,pcs_sorted_off));
hold on
line([32,32],[0,2000],'Color','w')
line([72,72],[0,2000],'Color','w')
line([0,2000],[32,32],'Color','w')
line([0,2000],[72,72],'Color','w')
ylabel('Odd')
xlabel('OFF');
title('Odd vs OFF');
yticks([1,32,72,105]);
yticklabels({'0','50','110','157'});
xticks([1,32,72,105]);
xticklabels({'0','50','110','157'});

set(gca,'FontSize',16)


%% Plot correlation but using all ON laps

[pcs_sorted_on, sort_ind] = sb.sort.rasterMapByPeak(pcs_avg_dff_on,'normalize',true);
pcs_sorted_off = sb.sort.rasterMapByPeak(pcs_avg_dff_off,'sort_order',sort_ind,'normalize',true);

% Normalize the cells to each other
%pcs_sorted_on = zScoreNormalize(pcs_sorted_on,'row');
%pcs_sorted_off = zScoreNormalize(pcs_sorted_off,'row');

figure(3);
clf;
subplot(2,2,1);
imagesc(pcs_sorted_on,[0,4]);
title('ON laps');
xticks([1,32,72,90,105]);
xticklabels({'0','50','110','140','157'})
ylabel('Cell #');
line([32,32],[0,2000],'Color','w')
line([72,72],[0,2000],'Color','w')
set(gca,'FontSize',16)

subplot(2,2,2);
imagesc(pcs_sorted_off,[0,4]);
title('OFF laps');
xticks([1,32,72,90,105]);
xticklabels({'0','50','110','140','157'})
ylabel('Cell #');
line([32,32],[0,2000],'Color','w')
line([72,72],[0,2000],'Color','w')
set(gca,'FontSize',16)

subplot(2,2,3);
imagesc(corr(pcs_sorted_on,pcs_sorted_off))
hold on;
line([32,32],[0,105])
line([72,72],[0,105])
title('Correlation heatmap');
set(gca,'FontSize',16)

% Compare the correlation by dividing the track into segments
correlation = diag(corr(pcs_sorted_on,pcs_sorted_off));
start_mean = nanmean(correlation(1:31));
start_sem = nanstd(correlation(1:31))*2;
first_mean = nanmean(correlation(32:71));
first_sem = nanstd(correlation(32:71))*2;
second_mean = nanmean(correlation(72:105));
second_sem = nanstd(correlation(72:105))*2;

subplot(2,2,4);
bar([start_mean,first_mean, second_mean])
hold on;
errorbar(1,start_mean,-start_sem,start_sem,'Color','k','LineWidth',1.5)
errorbar(2,first_mean,-first_sem,first_sem,'Color','k','LineWidth',1.5)
errorbar(3,second_mean,-second_sem,second_sem,'Color','k','LineWidth',1.5)

xticklabels({'Track start','First landmark','Second landmark','Track end'})
xtickangle(25)
ylim([0,1])
ylabel('Correlation');
title('Average correlation in track segments');
set(gca,'FontSize',16)

%% Plot the max correlation values along the position map

figure(2); clf;
plot(nanmax(corr(pcs_sorted_on_odd,pcs_sorted_on_even)),'Color','k');
hold on
plot(nanmax(corr(pcs_sorted_on_odd,pcs_sorted_off)),'Color','r');
ylim([0.5,1])


%% Plot the distribution of peak position for place cells in even vs odd vs OFF laps.

[max_amplitude_odd,max_position_odd] = max(pcs_avg_dff_on_odd');
[max_amplitude_even,max_position_even] = max(pcs_avg_dff_on_even');
[max_amplitude_off,max_position_off] = max(pcs_avg_dff_off');

figure(2); clf;
subplot(1,2,1)
histogram(max_position_odd,20);
hold on
histogram(max_position_off,20);

subplot(1,2,2);
boxplot([max_amplitude_odd',max_amplitude_even',max_amplitude_off']);


%% Plot the shift in place field max position

figure(3);clf; 
histogram(max_position_odd - max_position_even,210); 
hold on; 
histogram(max_position_odd - max_position_off,210)



%% Plot all place cells sorted for even, odd and OFF laps

% Normalize the cells to each other
pcs_sorted_on_odd = zScoreNormalize(pcs_sorted_on_odd,'row');
pcs_sorted_on_even = zScoreNormalize(pcs_sorted_on_even,'row');
pcs_sorted_off = zScoreNormalize(pcs_sorted_off,'row');

figure(3);
clf;
subplot(1,3,1);
imagesc(pcs_sorted_on_odd,[0,4]);
title('ON laps (Odd)');
xticks([1,32,72,90,105]);
xticklabels({'0','50','110','140','157'})
ylabel('Cell #');
line([32,32],[0,2000],'Color','w')
line([72,72],[0,2000],'Color','w')

set(gca,'FontSize',16)

subplot(1,3,2);
imagesc(pcs_sorted_on_even,[0,4]);
title('ON laps (Even)');
xticks([1,32,72,90,105]);
xticklabels({'0','50','110','140','157'})
ylabel('Cell #');
line([32,32],[0,2000],'Color','w')
line([72,72],[0,2000],'Color','w')
set(gca,'FontSize',16)

subplot(1,3,3);
imagesc(pcs_sorted_off,[0,4]);
title('OFF laps');
xticks([1,32,72,90,105]);
xticklabels({'0','50','110','140','157'})
ylabel('Cell #');
line([32,32],[0,2000],'Color','w')
line([72,72],[0,2000],'Color','w')
set(gca,'FontSize',16)


%% What percentage of place cells loose tuning when the landmark is off?

% Find all cells that have tuning around the first landmark
[~,position] = max(pcs_avg_dff_on_odd');
peak_positions =position(find([position>32] & [position<70]));
relevant_place_cells_even = pcs_avg_dff_on_even(find([position>32] & [position<70]),:);
relevant_place_cells_odd = pcs_avg_dff_on_odd(find([position>32] & [position<70]),:);
relevant_place_cells_off = pcs_avg_dff_off(find([position>32] & [position<70]),:);

on_response_odd = [];
on_response_even = [];
off_response = [];
window_size = 5;
for c = 1:length(peak_positions)
    
   on_response_even(c) = nanmean(relevant_place_cells_even(c,peak_positions(c)-window_size:peak_positions(c)+window_size-1));
   on_response_odd(c) = nanmean(relevant_place_cells_odd(c,peak_positions(c)-window_size:peak_positions(c)+window_size-1));
   off_response(c) = nanmean(relevant_place_cells_off(c,peak_positions(c)-window_size:peak_positions(c)+window_size-1));
    
end

% Plot histogram of the difference in response in on vs off laps, with even
% vs odd as baseline.
figure(1); clf;
histogram((on_response_even-on_response_odd)./(on_response_even+on_response_odd),20,'FaceColor','k','BinEdges',-1:0.1:1)
hold on;
histogram((off_response-on_response_odd)./(off_response+on_response_odd),20,'FaceColor','b','BinEdges',-1:0.1:1)
legend({'Even vs Odd','OFF vs O dd'})
ylabel('# of cells');
xlabel('% change in amplitude')
xticks([-1,-0.5,0,0.5,1]);
xticklabels({'-100','-50','0','50','100'})

set(gca,'FontSize',16)

% Plot the response amplitude change as a function of distance from the
% landmark that is manipulated
[distance_from_landmark_sorted,sort_ind] = sort(peak_positions-32);

on_response_even_sorted_by_distance = on_response_even(sort_ind);
on_response_odd_sorted_by_distance = on_response_odd(sort_ind);
off_response_sorted_by_distance = off_response(sort_ind);

% Compute the amplitude change for each roi in even vs odd and off vs odd
even_odd_responses = (on_response_even_sorted_by_distance-on_response_odd_sorted_by_distance)./(on_response_even_sorted_by_distance+on_response_odd_sorted_by_distance);
off_odd_responses = (off_response_sorted_by_distance-on_response_odd_sorted_by_distance)./(off_response_sorted_by_distance+on_response_odd_sorted_by_distance);

% Discretize the distance vector
distance_from_landmark_sorted = discretize(distance_from_landmark_sorted,11); 

value_off = [];
value_on = [];
for x = 1:max(distance_from_landmark_sorted)
   value_off(x) = nanmean(off_odd_responses(find([distance_from_landmark_sorted == x]))) ;
   value_on(x) = nanmean(even_odd_responses(find([distance_from_landmark_sorted == x]))) ;
end

figure(2);clf;
plot(value_on,'Color','k');
hold on
plot(value_off);

ylim([-1,0])



%% Place cells: Look at the correlation between each bin across on off laps
odd_even_corr = [];
odd_off_corr = [];

% % bin the data
% start_inds = 1:5:105;
% pcs_sorted_on_odd_binned = [];
% for i = 1:21
%     pcs_sorted_on_odd_binned(i,1) = pcs_sorted_on_odd(start_inds(i):start_inds(i)+4,1);   
% end

for i = 1:105
    temp = corrcoef(pcs_sorted_on_odd(:,i),pcs_sorted_on_even(:,i));
    odd_even_corr(i) = temp(1,2);
    temp = corrcoef(pcs_sorted_on(:,i),pcs_sorted_off(:,i));
    odd_off_corr(i) = temp(1,2);
end

figure(3); clf; 
plot(odd_even_corr);
hold on;
plot(odd_off_corr)
line([32,32],[0,1]);
line([72,72],[0,1]);
ylim([0,1])
title('Population correlation place cells')
set(gca,'FontSize',16)
xlim([0,105])



