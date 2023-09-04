
% The aim is to see if we can train a position classifier on the 'landmark
% ON (odd)' laps and see how well it performs on 'landmark OFF' laps vs
% 'landmark ON (even)' laps as a control.

% Position decoding is the same Bayesian method used in Mao et al (2020;
% Current Biology).

%% ----- Load data -----
% Clean up workspace
clearvars -except mData sData data

% Parameters
n_bins_to_pool = 15;
tau = n_bins_to_pool/31;
rng(123) % set random number generator seed

mean_decoding_errors_ON = [];
mean_decoding_errors_OFF = [];
decoding = struct();

for session = 1:10
  
    % Try first by only using one session
    sData = mData(session).sData;

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
    deconv = deconv + randi(10,size(deconv,1),size(deconv,2))*0.0000000000000000001;

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
        deconv_test_binned(:,i) = nanmean(deconv_test(:,s:s+n_bins_to_pool-1),2);
        pos_test_binned(:,i) = nanmean(pos_test(:,s:s+n_bins_to_pool-1),2);
        lap_test_binned(:,i) = nanmean(lap_test(:,s:s+n_bins_to_pool-1),2);
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

    mean_decoding_errors_ON(session) = mean_decoding_error_val;
    mean_decoding_errors_OFF(session) = mean_decoding_error_test;
    
    decoding(session).decoding_error_ON = decoding_error_val*1.5;
    decoding(session).decoding_error_OFF = decoding_error_test*1.5;
    
end


disp('done')

%
figure(2);clf;
boxplot([mean_decoding_errors_ON;mean_decoding_errors_OFF]')
ylim([0,35])
ylabel('Mean decoding error (cm)');
xticklabels({'Landmark ON','Landmark OFF'})
title('Decoding error')
set(gca,'FontSize',16)


%
% Plot result
figure(1); clf;
subplot(2,1,1);
plot(pos_val_binned,'Color','k','LineWidth',1.5);
hold on;
plot(predicted_pos_val,'Color','r','LineWidth',1.5);
set(gca,'FontSize',16)
text(10,110,sprintf('Mean decoding error: %.1f cm',mean_decoding_error_val),'FontSize',12)
xlim([0,length(predicted_pos_val)])
title(sprintf('Validation - %s',sData.sessionInfo.sessionID(1:17)))

xticks([1,length(predicted_pos_val)])
%xticklabels({0,sprintf('%.1f',(length(predicted_pos_val)*(n_bins_to_pool/31))/60)})
xlabel('Time (min)');
ylabel('Position (cm)');

legend({'Actual position','Predicted position'})

subplot(2,1,2);
plot(pos_test_binned,'Color','k','LineWidth',1.5);
hold on;
plot(predicted_pos_test,'Color','r','LineWidth',1.5);
set(gca,'FontSize',16)
text(10,110,sprintf('Mean decoding error: %.1f cm',mean_decoding_error_test),'FontSize',12)
xlim([0,length(predicted_pos_test)])
title(sprintf('Testing - %s',sData.sessionInfo.sessionID(1:17)))

xticks([1,length(predicted_pos_test)])
%xticklabels({0,sprintf('%.1f',(length(predicted_pos_test)*2)/60)})
xlabel('Time (min)');
ylabel('Position (cm)');

legend({'Actual position','Predicted position'})


%% Plot the decoding errors as a function of time
figure(1);clf;

subplot(2,1,1);
plot(decoding_error_val)
hold on
line([0,10000],[nanmean(decoding_error_val),nanmean(decoding_error_val)],'LineWidth',2)
xlim([0,length(decoding_error_val)])
title('Validation');
xticks([1,length(decoding_error_val)])
ylabel('Error (cm)');
set(gca,'FontSize',16)

subplot(2,1,2);
plot(decoding_error_test)
hold on
line([0,10000],[nanmean(decoding_error_test),nanmean(decoding_error_test)],'LineWidth',2)
xlim([0,length(decoding_error_test)])
title('Testing');
xticks([1,length(decoding_error_test)])
xlabel('Time (bins)');
ylabel('Error (cm)');
set(gca,'FontSize',16)


%% Plot the decoding errors as a function of position on track
[~,~,unique_laps] = unique(round(lap_val_binned));
pos = discretize(round(pos_val_binned),105);
rmap = sb.generate.rasterMaps(pos,unique_laps',decoding_error_val);
figure(1); clf;
plot(nanmean(rmap));


%% Compare the decoding error across all sessions for the area around the landmark which is omitted
decoding_error_val_at_landmark = [];
for i = 1:length(pos_val_binned)
    if ismember(round(pos_val_binned(i)),32:52)
        decoding_error_val_at_landmark = [decoding_error_val_at_landmark, decoding_error_val(i)];
    end
end

decoding_error_test_at_landmark = [];
for i = 1:length(pos_test_binned)
    if ismember(round(pos_test_binned(i)),32:52)
        decoding_error_test_at_landmark = [decoding_error_test_at_landmark, decoding_error_test(i)];
    end
end

figure(3); clf;
nanmean(decoding_error_test_at_landmark)
nanstd(decoding_error_test_at_landmark)
nanmean(decoding_error_val_at_landmark)

%% Show decoding for all sessions and the distributions
figure(4); clf;
subplot(1,4,1:3);
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

    plot([s,s],[nanmean(decoding(s).decoding_error_OFF)],'LineStyle','none','Marker','.','MarkerSize',20,'Color','r')
    hold on
    errorbar(s,decoding_mean,-decoding_sem,decoding_sem,'Color','r')

    
end

ylim([0,40])
xlim([0,11])
xticks([1:10])
xlabel('Session #');
set(gca,'FontSize',16)
ylabel('Decoding error (cm)');
subplot(1,4,4);

decoding_mean = nanmean(decoding_mean_ON_all);
decoding_sem = nanstd(decoding_mean_ON_all)/sqrt(length(decoding_mean_ON_all));
bar(1,nanmean(decoding_mean_ON_all),'FaceColor','k');
hold on
errorbar(1,decoding_mean,-decoding_sem,decoding_sem,'Color',[0.5,0.5,0.5],'LineWidth',2);

decoding_mean = nanmean(decoding_mean_OFF_all);
decoding_sem = nanstd(decoding_mean_OFF_all)/sqrt(length(decoding_mean_OFF_all));
bar(2,nanmean(decoding_mean_OFF_all),'FaceColor','r');
errorbar(2,decoding_mean,-decoding_sem,decoding_sem,'Color',[0.5,0.5,0.5],'LineWidth',2);
xticks([1,2])
xticklabels({'Control','OFF'})
xtickangle(25)
ylim([0,40])
set(gca,'FontSize',16)



%%
figure(2);clf;
plot([1.0025,1.0025],[2],'LineStyle','none','Marker','.','MarkerSize',20,'Color','k')
hold on
errorbar(1,2,-0.2,0.2,'Color','k')
ylim([0,30])
xlim([0,2])


