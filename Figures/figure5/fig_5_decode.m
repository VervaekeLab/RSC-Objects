%% ----- Load data -----
% Clean up workspace
clearvars -except mData sData data

% Parameters
n_bins_to_pool = 10;
tau = n_bins_to_pool/31;
rng(123) % set random number generator seed

decoding = [];

for session = 1:10
    % Try first by only using one session
    sData = mData(session).sData;
    if session == 1
        decoding = decodeSession(sData);
    else
        decoding(session) = decodeSession(sData);
    end
end

disp('done')

% Plot mean decoding error across all sessions
mean_decoding_errors_ON = [];
mean_decoding_errors_OFF = [];
for s = 1:length(decoding)
    mean_decoding_errors_ON = [mean_decoding_errors_ON,nanmean(decoding(s).decoding_error_ON)];
    mean_decoding_errors_OFF = [mean_decoding_errors_OFF,nanmean(decoding(s).decoding_error_OFF)];
end

% figure(2);clf;
% boxplot([mean_decoding_errors_ON;mean_decoding_errors_OFF]')
% ylim([0,35])
% ylabel('Mean decoding error (cm)');
% xticklabels({'Landmark ON','Landmark OFF'})
% title('Decoding error')
% set(gca,'FontSize',16)

% Plot example of decoding
% Showing 3 minute of the recording, which is 124 bins when using ~0.5s bin width
s = 6;
first_bin = 1;
last_bin = 258;
n_minutes = 2;
figure(1); clf;
subplot(2,1,1);
line_width = 0.5;
plot(decoding(s).pos_val_binned(first_bin:last_bin),'Color','k','LineWidth',line_width);
hold on;
plot(decoding(s).predicted_pos_val(first_bin:last_bin),'Color','r','LineWidth',line_width);
set(gca,'FontSize',16)
text(10,110,sprintf('Mean decoding error: %.1f cm',nanmean(decoding(s).decoding_error_ON)),'FontSize',12)
xlim([0,length(decoding(s).predicted_pos_val(first_bin:last_bin))])
title('Control')

xticks([1,length(decoding(s).predicted_pos_val(first_bin:last_bin))])
xticklabels({0,n_minutes})
ylim([-4,109])
yticks([1,105])
yticklabels({0,157})
xlabel('Time (min)');
ylabel('Position (cm)');

legend({'Actual position','Predicted position'})

subplot(2,1,2);
plot(decoding(s).pos_test_binned(first_bin:last_bin),'Color','k','LineWidth',line_width);
hold on;
plot(decoding(s).predicted_pos_test(first_bin:last_bin),'Color','r','LineWidth',line_width);
set(gca,'FontSize',16)
text(10,110,sprintf('Mean decoding error: %.1f cm',nanmean(decoding(s).decoding_error_OFF)),'FontSize',12)
xlim([0,length(decoding(s).predicted_pos_test(first_bin:last_bin))])
title('Landmark omitted')
xticks([1,length(decoding(s).predicted_pos_test(first_bin:last_bin))])
xticklabels({0,n_minutes})
ylim([-4,109])
yticks([1,105])
yticklabels({0,157})


xlabel('Time (min)');
ylabel('Position (cm)');

%legend({'Actual position','Predicted position'})

set(gcf,'renderer', 'painters', 'Position', [2000 200 700 300])

%% Show decoding for all sessions and the distributions
figure(3); clf;
subplot(1,4,1:2);
i = 1;
decoding_mean_ON_all = [];
decoding_mean_OFF_all = [];
for s = 1:10
    
    % Plot the decode ON
    decoding_mean = [nanmean(decoding(s).decoding_error_ON)];
    decoding_mean_ON_all = [decoding_mean_ON_all,decoding_mean];
    decoding_sem = nanmean(decoding(s).decoding_error_ON)/sqrt(length(decoding(s).decoding_error_ON));
    plot([s,s],decoding_mean,'LineStyle','none','Marker','.','MarkerSize',20,'Color','k')
    hold on
    errorbar(s,decoding_mean,-decoding_sem,decoding_sem,'Color','k')
    
    
    % Plot the decode ON
    decoding_mean = [nanmean(decoding(s).decoding_error_OFF)];
    decoding_mean_OFF_all = [decoding_mean_OFF_all,decoding_mean];
    decoding_sem = nanmean(decoding(s).decoding_error_OFF)/sqrt(length(decoding(s).decoding_error_OFF));

    plot([s,s],[nanmean(decoding(s).decoding_error_OFF)],'LineStyle','none','Marker','.','MarkerSize',20,'Color','b')
    hold on
    errorbar(s,decoding_mean,-decoding_sem,decoding_sem,'Color','b')
    
end

ylim([0,40])
xlim([0,11])
xticks([1:10])
xlabel('Session #');
set(gca,'FontSize',16)
ylabel('Decoding error (cm)');
subplot(1,4,3:4);

decoding_mean = nanmean(decoding_mean_ON_all);
decoding_sem = nanstd(decoding_mean_ON_all)/sqrt(length(decoding_mean_ON_all));
bar(1,nanmean(decoding_mean_ON_all),'FaceColor','k');
hold on
errorbar(1,decoding_mean,-decoding_sem,decoding_sem,'Color',[0.5,0.5,0.5],'LineWidth',2);

decoding_mean = nanmean(decoding_mean_OFF_all);
decoding_sem = nanstd(decoding_mean_OFF_all)/sqrt(length(decoding_mean_OFF_all));
bar(2,nanmean(decoding_mean_OFF_all),'FaceColor','b');
errorbar(2,decoding_mean,-decoding_sem,decoding_sem,'Color',[0.5,0.5,0.5],'LineWidth',2);
xticks([1,2])
xticklabels({'Control','Omitted'})
xtickangle(25)
ylim([0,40])
set(gca,'FontSize',16)
set(gcf,'renderer', 'painters', 'Position', [2000 200 500 300])









%% ---- HELPER FUNCTION ----- 
function decoding_results = decodeSession(sData)

    % Parameters
    n_bins_to_pool = 6;
    tau = n_bins_to_pool/31;
    rng(123) % set random number generator seed

    mean_decoding_errors_ON = [];
    mean_decoding_errors_OFF = [];
    
    % init    
    decoding_results = struct();
    
    % Extract data
    deconv = double(sData.imdata.roiSignals(2).deconv(:,:));
    pos = sData.behavior.wheelPosDsBinned;
    lap = sData.behavior.wheelLapDsBinned;
    speed = sData.behavior.runSpeedDs;

    % Exclude silent cells
    deconv = deconv(find(sum(deconv,2)),:);

    % Remove all bins where running speed is < 1 cm.
    bins_to_include = find([speed<1] == 0);
    deconv = deconv(:,bins_to_include);
    pos = pos(bins_to_include);
    lap = lap(bins_to_include);
    speed = speed(bins_to_include);

    % Add non zero noise to deconv
    deconv = deconv + randi(10,size(deconv,1),size(deconv,2))*0.0000000000000001;

    % Smooth the deconv signal (Mao does this in his code, but don't write so in the paper) 
    %deconv = smoothdata(deconv,2,'gaussian',10);

    % Shift lap to start on 1, not -1 or 0.
    lap = lap + (abs(min(lap))*2);

    % Split into train and test set (odd vs even laps)
    laps_landmark_on = find(sData.landmarks.motor.on(:,1) == 1);
    laps_landmark_off = find(sData.landmarks.motor.on(:,1) == 0);

    laps_landmark_on_even = laps_landmark_on(find(mod(laps_landmark_on,2) == 0));
    laps_landmark_on_odd = laps_landmark_on(find(mod(laps_landmark_on,2) == 1));

    % Find indices for the OFF laps, ON (even) and ON (odd) laps
    off_laps = ismember(lap,laps_landmark_off);
    on_odd_laps = ismember(lap,laps_landmark_on_odd);
    on_even_laps = ismember(lap,laps_landmark_on_even);

    % Divide the data into train, val, test
    deconv_train = deconv(:,on_odd_laps);
    deconv_val = deconv(:,on_even_laps);
    deconv_test = deconv(:,off_laps);
    pos_train = pos(on_odd_laps);
    pos_test = pos(off_laps);
    pos_val = pos(on_even_laps);
    lap_train = lap(on_odd_laps);
    lap_test = lap(off_laps);
    lap_val = lap(on_even_laps);

    % ----- Compute position tuning curves (1D) -----
    % Compute trial averaged and occupancy-normalized position tuning curves
    % for all neurons, using train dataset

    % Occupancy 
    occupancy = zeros(1,105);
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
        cell_tuning = zeros(1,105);

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
    lap_test_binned = [];
    lap_val_binned = [];

    % Test data (OFF laps)
    i = 1;
    for s = 1:n_bins_to_pool:length(deconv_test)-n_bins_to_pool
%         deconv_test_binned(:,i) = nanmean(deconv_test(:,s:s+n_bins_to_pool-1),2);
%         pos_test_binned(:,i) = nanmean(pos_test(:,s:s+n_bins_to_pool-1),2);
%         lap_test_binned(:,i) = nanmean(lap_test(:,s:s+n_bins_to_pool-1),2);
        deconv_test_binned(:,i) = nanmean(deconv_test(:,s:s+n_bins_to_pool-1),2);
        pos_test_binned(:,i) = pos_test(:,s+n_bins_to_pool-1);
        lap_test_binned(:,i) = lap_test(:,s+n_bins_to_pool-1);
        i = i+1;
    end

    % Validation data (ON even laps)
    i = 1;
    for s = 1:n_bins_to_pool:length(deconv_val)-n_bins_to_pool
        deconv_val_binned(:,i) = nanmean(deconv_val(:,s:s+n_bins_to_pool-1),2);
        pos_val_binned(:,i) = nanmean(pos_val(:,s:s+n_bins_to_pool-1),2);
        lap_val_binned(:,i) = nanmean(lap_val(:,s:s+n_bins_to_pool-1),2);
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
    decoding_error_val = abs(pos_val_binned-predicted_pos_val);

    % Correct for circularity
    need_correcting = find(decoding_error_test>(105/2));
    decoding_error_test(need_correcting) = 105-decoding_error_test(need_correcting);
    mean_decoding_error_test = nanmean(decoding_error_test)*1.5;

    need_correcting = find(decoding_error_val>(105/2));
    decoding_error_val(need_correcting) = 105-decoding_error_val(need_correcting);
    mean_decoding_error_val = nanmean(decoding_error_val)*1.5;
    
    decoding_results.decoding_error_ON = decoding_error_val*1.5;
    decoding_results.decoding_error_OFF = decoding_error_test*1.5;
   
    decoding_results.pos_val_binned = pos_val_binned;
    decoding_results.predicted_pos_val = predicted_pos_val;
    
    decoding_results.pos_test_binned = pos_test_binned;
    decoding_results.predicted_pos_test = predicted_pos_test;
    
    decoding_results.pos_test_original = pos_test;

end
