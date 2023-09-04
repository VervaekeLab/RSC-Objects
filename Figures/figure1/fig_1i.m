%% Control experiments to see how many "landmark cells" is present when mice have had their whiskers cut

disp("Loading data ...")

% Load main dataset
[mData,data] = sb.database.load.rsc_manipulation;

% Load control dataset
[mData2,data2] = sb.database.load.rsc_whisker_control;

% Load blank wheel dataset
[mData3,data3] = sb.database.load.rsc_blank_control;

disp("All data is loaded.")

%% Compare the percentage of landmark cells between main dataset and control
landmark_cell_color = [239,35,60]/255;
place_cell_color = [51,51,51]/255;

% Get the percentage of cells
main_lcs_pct = (data.n_lcs + data.n_tcs) ./ data.n_acs;
control_lcs_pct = data2.n_lcs ./ data2.n_acs;
blank_lcs_pct = data3.n_lcs./ data3.n_acs;

main_lcs_stderror = std(main_lcs_pct) / sqrt(length(main_lcs_pct));
control_lcs_stderror = std(control_lcs_pct) / sqrt(length(control_lcs_pct));
blank_lcs_stderror = std(blank_lcs_pct) / sqrt(length(blank_lcs_pct));

% Plot
figure(1);
clf;
subplot(1,2,2)

% Main data
main_lcs_avg = nanmean(main_lcs_pct);
low_error = -main_lcs_stderror;
high_error = main_lcs_stderror;

b = bar(1,main_lcs_avg);
b.FaceColor = landmark_cell_color;
hold on
er = errorbar(1,main_lcs_avg,low_error, high_error);
er.Color = [0,0,0];
er.LineStyle = 'none';
er.LineWidth = 2;

% Control data
control_lcs_avg = nanmean(control_lcs_pct);
low_error = -control_lcs_stderror;
high_error = control_lcs_stderror;

b = bar(2,control_lcs_avg);
b.FaceColor = [0.8,0.8,0.8];
hold on
er = errorbar(2,control_lcs_avg,low_error, high_error);
er.Color = [0,0,0];
er.LineWidth = 2;
er.LineStyle = 'none';
ylim([0,0.15])

% Blank data
blank_lcs_avg = nanmean(blank_lcs_pct);
low_error = -blank_lcs_stderror;
high_error = blank_lcs_stderror;

b = bar(3,blank_lcs_avg);
b.FaceColor = [0.8,0.8,0.8];
hold on
er = errorbar(3,blank_lcs_avg,low_error, high_error);
er.Color = [0,0,0];
er.LineWidth = 2;
er.LineStyle = 'none';
ylim([0,0.15])

yticks([0,0.05,0.1,0.15])
yticklabels({'0','5','10','15'})
ylabel('% of cells');
xticks([1,2,3])
xticklabels({'Whiskers','No whiskers','Blank'})
xtickangle(25)
title('Landmark cells');
set(gca,'FontSize',18);


% Compare the percentage of place cells between main dataset and control

% Get the percentage of cells
main_pcs_pct = data.n_pcs ./ data.n_acs;
control_pcs_pct = data2.n_pcs ./ data2.n_acs;
blank_pcs_pct = data3.n_pcs ./ data3.n_acs;

main_pcs_stderror = std(main_pcs_pct) / sqrt(length(main_pcs_pct));
control_pcs_stderror = std(control_pcs_pct) / sqrt(length(control_pcs_pct));
blank_pcs_stderror = std(blank_pcs_pct) / sqrt(length(blank_pcs_pct));

% Plot
figure(1);
subplot(1,2,1)

% Main data
main_pcs_avg = nanmean(main_pcs_pct);
low_error = -main_pcs_stderror;
high_error = main_pcs_stderror;

b = bar(1,main_pcs_avg);
b.FaceColor = place_cell_color;
hold on
er = errorbar(1,main_pcs_avg,low_error, high_error);
er.Color = [0,0,0];
er.LineStyle = 'none';
er.LineWidth = 2;
ylim([0,0.12])

% Control data
control_pcs_avg = nanmean(control_pcs_pct);
low_error = -control_pcs_stderror;
high_error = control_pcs_stderror;

b = bar(2,control_pcs_avg);
b.FaceColor = [0.8,0.8,0.8];
hold on
er = errorbar(2,control_pcs_avg,low_error, high_error);
er.Color = [0,0,0];
er.LineWidth = 2;
er.LineStyle = 'none';
ylim([0,0.4])

% Blank data
blank_pcs_avg = nanmean(blank_pcs_pct);
low_error = -control_pcs_stderror;
high_error = control_pcs_stderror;

b = bar(3,blank_pcs_avg);
b.FaceColor = [0.8,0.8,0.8];
hold on
er = errorbar(3,blank_pcs_avg,low_error, high_error);
er.Color = [0,0,0];
er.LineWidth = 2;
er.LineStyle = 'none';
ylim([0,0.4])

yticks([0,0.1,0.2,0.3,0.4])
yticklabels({'0','10','20','30','40'})
ylabel('% of cells');
xticks([1,2,3])
xticklabels({'Whiskers','No whiskers','Blank'})
xtickangle(25)
title('Place cells');
set(gca,'FontSize',18);

% Set correct size of the figure
set(gcf,'renderer', 'painters', 'Position', [1500,100,800,350])

%% Visualize all landmark cells that is found
avg_tuning_lcs = [];
for s = 1:length(mData2)
    
    try
        if isempty(avg_tuning_lcs)
            avg_tuning_lcs = sb.generate.averagedTuningFromRasterMaps(mData2(s).rmaps.lcs.dff_motor_on);
        else
            avg_tuning_lcs = [avg_tuning_lcs; sb.generate.averagedTuningFromRasterMaps(mData2(s).rmaps.lcs.dff_motor_on)];
        end
    end
    
    
end

sb.plot.rasterMaps(sb.sort.rasterMapByPeak(avg_tuning_lcs));

%% Visualize all place cells
avg_tuning_pcs = [];
for s = 1:length(mData2)
    
    try
        if isempty(avg_tuning_pcs)
            avg_tuning_pcs = sb.generate.averagedTuningFromRasterMaps(mData2(s).rmaps.pcs.dff_motor_on);
        else
            avg_tuning_pcs = [avg_tuning_pcs; sb.generate.averagedTuningFromRasterMaps(mData2(s).rmaps.pcs.dff_motor_on)];
        end
    end
    
end

sb.plot.rasterMaps(sb.sort.rasterMapByPeak(avg_tuning_pcs));


%% Create bar plot of percentage of cells
pct_lcs = [];
for s = 1:length(mData2)
    pct_lcs(s) = length(mData2(s).rmaps.lcs.ind)/size(mData2(s).rmaps.dff,3);
end

figure(10);clf;
bar(pct_lcs);
ylim([0,0.1])


