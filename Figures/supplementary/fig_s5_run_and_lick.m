%% Load data
clear; 
[mData, data] = sb.database.load.rsc_manipulation;

%% Create figure
line_width = 2; 
% Get the rastermap of running speed for each session
for s = 1:length(mData)    
    
    % Extract timeseries data
    pos = mData(s).sData.behavior.wheelPosDsBinned;
    lap = mData(s).sData.behavior.wheelLapDsBinned;
    run = mData(s).sData.behavior.runSpeedDs;
    lick = mData(s).sData.daqdata.lickSignal(mData(s).sData.daqdata.frameIndex);
    
    % Create rastermap of the running speed binned along the track
    running_rmap = sb.generate.rasterMaps(pos,lap,run,'smoothing',1);
    lick_rmap = sb.generate.rasterMaps(pos,lap,lick,'smoothing',1);
    
    % Find indices for laps when the landmark is on or off.
    on_laps = find([mData(s).sData.landmarks.motor.on(:,1) == 1]);
    off_laps = find([mData(s).sData.landmarks.motor.on(:,1) == 0]);
    
    if s == 1
        running_rmap_landmark_on = running_rmap(on_laps,:);
        running_rmap_landmark_off = running_rmap(off_laps,:);
        
        lick_rmap_landmark_on = lick_rmap(on_laps,:);
        lick_rmap_landmark_off = lick_rmap(off_laps,:);
        
        
    else
        running_rmap_landmark_on = cat(1,running_rmap_landmark_on,running_rmap(on_laps,:));
        running_rmap_landmark_off = cat(1,running_rmap_landmark_off,running_rmap(off_laps,:));
        
        lick_rmap_landmark_on = cat(1,lick_rmap_landmark_on,lick_rmap(on_laps,:));
        lick_rmap_landmark_off = cat(1,lick_rmap_landmark_off,lick_rmap(off_laps,:));
    end
    
   
end

% Plot individual runs and the average
figure(1); clf;
color_map = flipud(cbrewer('div','PuOr',5));
color_map(5,:) = [0,0,0]; % Set the off color to black
color_alpha = 0.1;
line_width = 1;
subplot(2,1,1);

% Plot speed when landmark is present
ops = struct();
ops.error = 'std';
ops.handle     = figure(1);
ops.alpha      = color_alpha;
ops.line_width = line_width;
ops.color_area = color_map(1,:);
ops.color_line = color_map(1,:);
plot_areaerrorbar(running_rmap_landmark_on,ops);

hold on;

% Plot speed when landmark is omitted
ops = struct();
ops.error = 'std';
ops.handle     = figure(1);
ops.alpha      = color_alpha;
ops.line_width = line_width;
ops.color_area = color_map(5,:);
ops.color_line = color_map(5,:);
plot_areaerrorbar(running_rmap_landmark_off,ops);

% Draw lines indicating the landmarks
line([32,32],[-100,100],'Color','r')
line([72,72],[-100,100],'Color','r')
line([92,92],[-100,100],'Color','k')

%legend({'Landmark present','Landmark omitted'})
xticks([0,105])
xlim([0,105])
xticklabels({0,157})
title('Running speed');
ylabel('Running speed (cm/s)');
xlabel('Position (cm)');
ylim([0,73])
yticks([0,70])
set(gca,'FontSize',16);
box('off')

subplot(2,1,2);
% ops = struct();
% ops.error = 'var';
% ops.handle     = figure(1);
% ops.alpha      = 0.25;
% ops.line_width = 2;
% ops.color_area = [0.2,0.2,0.2];
% ops.color_line = [0,0,0];
% plot_areaerrorbar(lick_rmap_landmark_on,ops);


%plot(nanmean(lick_rmap_landmark_on,1),'Color','k','LineWidth',line_width);
%hold on
%plot(nanmean(lick_rmap_landmark_off,1),'Color','b','LineWidth',line_width);


lick_probability_on = nansum(lick_rmap_landmark_on,1);
lick_probability_on = lick_probability_on/sum(lick_probability_on);
lick_probability_off = nansum(lick_rmap_landmark_off,1);
lick_probability_off = lick_probability_off/sum(lick_probability_off);

area(nansum(lick_probability_on,1),'FaceColor',color_map(1,:),'FaceAlpha',color_alpha,'EdgeColor',color_map(1,:),'LineWidth',line_width);
hold on
area(nansum(lick_probability_off,1),'FaceColor',color_map(5,:),'FaceAlpha',color_alpha,'EdgeColor',color_map(5,:),'LineWidth',line_width);

% Add lines indicating the landmarks
line([32,32],[-100,100],'Color','r')
line([72,72],[-100,100],'Color','r')
line([92,92],[-100,100],'Color','k')

legend({'Landmark present','Landmark omitted'})

xticks([0,105])
xlim([0,105])
xticklabels({0,157})
title('Lick probability');
ylabel('Licks (A.U.)');
xlabel('Position (cm)');
ylim([0,0.07])
yticks([0,0.06])
set(gca,'FontSize',16);
box('off')

% Set figure size
set(gcf,'renderer', 'painters', 'Position', [100 200 500 700])
