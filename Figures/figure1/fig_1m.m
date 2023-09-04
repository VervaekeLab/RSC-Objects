
%% Compute landmark preference ratio 
amplitude_change = [];
amplitude_change_control = [];
amplitude_change_low_speed = [];
amplitude_change_high_speed = [];
avg_tuning_map_on = [];
avg_tuning_map_off = [];
amplitude_session = [];
amplitude_roi_index = [];
pval = [];
for s = 1:length(mData)
    
    % Include only landmark cells
    cells_to_include = [mData(s).rmaps.tcs.ind,mData(s).rmaps.lcs.ind]; 
    
    % Create averge tuning maps
    avg_tuning_map_on = cat(1,avg_tuning_map_on,sb.generate.averagedTuningFromRasterMaps(mData(s).rmaps.deconv_motor_on(:,:,cells_to_include)));
    avg_tuning_map_off = cat(1,avg_tuning_map_off,sb.generate.averagedTuningFromRasterMaps(mData(s).rmaps.deconv_motor_off(:,:,cells_to_include)));
    
    % Compute the amplitude ratio between first and second landmark and as
    % a control
    amplitude_change = [amplitude_change, sb.quantify.landmarkPreferenceRatio(mData(s).rmaps.deconv_motor_on(:,:,cells_to_include))];
%     amplitude_change_control = [amplitude_change_control, sb.quantify.landmarkPreferenceRatioControl(mData(s).rmaps.deconv_motor_on(:,:,cells_to_include))];
    amplitude_shuffle = sb.quantify.landmarkPreferenceRatioShuffle(mData(s).rmaps.deconv_motor_on(:,:,cells_to_include));
    amplitude_change_control = [amplitude_change_control,amplitude_shuffle(:)'];

    % Add session number and cell indices of each that is included
    amplitude_session = [amplitude_session,ones(1,length(cells_to_include))*s];
    amplitude_roi_index = [amplitude_roi_index,cells_to_include];
    
    %%%%%%%% Not sure how useful these high/low speed stratified ratios are
    % Low speed laps
    laps_to_include = max(mData(s).rmaps.run_motor_on')<50;
    
    if sum(laps_to_include)>1
        amplitude_change_low_speed = [amplitude_change_low_speed,sb.quantify.landmarkPreferenceRatio(mData(s).rmaps.deconv_motor_on(laps_to_include,:,cells_to_include))];
    end
    
    % High speed laps
    laps_to_include = max(mData(s).rmaps.run_motor_on')>50; 
    if sum(laps_to_include)>1
        amplitude_change_high_speed = [amplitude_change_high_speed,sb.quantify.landmarkPreferenceRatio(mData(s).rmaps.deconv_motor_on(laps_to_include,:,cells_to_include))];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
end

[amplitude_change_sorted,sorted_amplitude_change_index] = sort(amplitude_change);

fprintf('Out of the %i landmark cells %.1f %% showed a biased greater than 100%% \n', length(amplitude_change),(sum(abs(amplitude_change)>100)/length(amplitude_change)*100))

disp("Finished loading data.")

amplitude_change(~isfinite(amplitude_change)) = [];
amplitude_change_control(~isfinite(amplitude_change_control)) = [];

figure(5); clf;

histogram(amplitude_change,'BinEdges',-350:10:350,'Normalization','probability','FaceColor','b');
hold on
ylim([0,0.1])
yticks([0,0.05,0.1])
yticklabels([0,5,10]);
ylabel('% of cells');
xlabel('% difference in amplitude');
xticks(-300:100:300)
text(0.1,0.8,'<--- Prefers first landmark','FontSize',15,'Units','normalized');
text(0.65,0.8,'Prefers second landmark --->','FontSize',15,'Units','normalized');
title('Amplitude difference between largest and smallest peak');
set(gca,'FontSize',16);

histogram(amplitude_change_control,'BinEdges',-350:10:350,'Normalization','probability','FaceColor','k','FaceAlpha',0.15);

set(gcf,'renderer', 'painters', 'Position', [2000,100,900,400])
% create line at 75th percentile of control group: xline(-79.941)
xline(quantile(amplitude_change_control,0.05))
xline(quantile(amplitude_change_control,0.95))
%% Show histogram of amplitude difference but split it into which is biggest, first or second landmark

figure(5); clf;
subplot(1,4,1:3);

%histogram(amplitude_change(amplitude_change>0),'BinEdges',0:10:350,'Normalization','probability','FaceColor','r');
%hold on;
%histogram(amplitude_change(amplitude_change<0),'BinEdges',-350:10:0,'Normalization','probability','FaceColor','b');

histogram(amplitude_change,'BinEdges',-350:10:350,'Normalization','probability','FaceColor','b');

ylim([0,0.08])
yticks([0,0.05])
yticklabels([0,5]);
ylabel('% of cells');
xlabel('% difference in amplitude');
xticks(-300:100:300)
text(-260,0.1,'<--- Prefers first landmark','FontSize',15);
text(50,0.1,'Prefers second landmark --->','FontSize',15);
title('Amplitude difference between largest and smallest peak');
set(gca,'FontSize',16);

% Show the box plot
subplot(1,4,4);
boxplot([abs(amplitude_change_control);abs(amplitude_change)]')
ylim([-20,500])
xticklabels({'Control','Difference'})
ylabel('% difference in amplitude')
title('Distribution compared with control');
set(gca,'FontSize',16)


% Set correct size of the figure
set(gcf,'renderer', 'painters', 'Position', [2000,100,1500,400])



%% Print some stats
% fprintf('First landmark: %.1f %% of cells have more than double the amplitude of second landmark\n', (sum(first_landmark_biased_SMI>100)/length(first_landmark_biased_SMI))*100)
% fprintf('Second landmark: %.1f %% of cells have more than double the amplitude of first landmark\n', (sum(second_landmark_biased_SMI>100)/length(second_landmark_biased_SMI))*100)

fprintf('In total %.1f %% of cells have more than twice the amplitude for one of the landmarks\n',(sum(abs(amplitude_change) > 100) / length(amplitude_change)) * 100)
fprintf('In total %.1f %% of cells have more than three times the amplitude for one of the landmarks\n',(sum(abs(amplitude_change) > 200) / length(amplitude_change)) * 100)

fprintf('In total %.1f %% of cells are below the 5th quantile of the shuffled distribution\n',(sum(amplitude_change < quantile(amplitude_change_control,0.05)) / length(amplitude_change)) * 100)
fprintf('In total %.1f %% of cells are above the 95th quantile of the shuffled distribution\n',(sum(abs(amplitude_change) > quantile(amplitude_change_control,0.95)) / length(amplitude_change)) * 100)


%% Is the difference significantly different to the control?

% We want to find out if the landmark cells consistently have a larger
% amplitude difference when comparing the first and second landmark versus
% a control. In this scenario, we choose the control to be odd vs even laps
% at the first landmark. We get the absolute value of difference in
% amplitude of both arrays and run a Wilcoxon sign rank test to test for
% consistency.

p_value = signrank(abs(amplitude_change),abs(amplitude_change_control));

prct_aboveThreshol = 100*length(find(abs(amplitude_change)> 79.941))./length(amplitude_change_control);


