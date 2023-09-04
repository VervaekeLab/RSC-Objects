%% Load the dataset for pupil sessions
clear;
[mData,data] = sb.database.load.rsc_pupil_sessions_50fps;

%% Run a rank sum test on the pupil radius events
ranksum_p = [];
rank_counter = 1;

for s = 1:3
    
    % Create rastermap for radius
    rmap_radius = sb.generate.rasterMaps(mData(s).sData.behavior.wheelPosDsBinned,mData(s).sData.behavior.wheelLapDsBinned,mData(s).sData.pupil.radius');
    rmap_radius = rmap_radius(1:50,:); % use only the dark laps
    
    % Detect all laps where there is a blink event on at first landmark
    events = [];
    counter = 1;
    for l = 1:size(rmap_radius,1)        
        baseline = nanmean(rmap_radius(l,11:31)) - (nanstd(rmap_radius(l,11:31))*2);
        
        if nanmin(rmap_radius(l,32:38))<baseline
            events(counter) = 1;
        else
            events(counter) = 0;
        end
        counter = counter + 1;
    end

    % Detect all laps where there is a blink event on at second landmark
    for l = 1:size(rmap_radius,1)        
        baseline = nanmean(rmap_radius(l,51:71)) - (nanstd(rmap_radius(l,51:71))*2);
        
        if nanmin(rmap_radius(l,72:78))<baseline
            events(counter) = 1;
        else
            events(counter) = 0;
        end

        counter = counter + 1;
    end

    % For each cell, get the neural response at landmark position 1
    for c = 1:length(mData(s).rmaps.level1.lcs)

        % Init emtpy array to contain the neural response
        nr = [];
        
        % Obtain rastermap for the current landmakr cell
        rmap = mData(s).rmaps.deconv(1:50,:,mData(s).rmaps.level1.lcs(c));
        
        % Get the average response at the landmark 
        counter = 1;
        for i = 1:size(rmap,1) % the first 50 laps
            nr(counter) = nanmax(rmap(i,32:42));
            counter = counter + 1;
            nr(counter) = nanmax(rmap(i,72:82));
            counter = counter + 1;
        end

        % Split nr responses into blink events
        nr_noevent = nr([events == 0]);
        nr_event = nr([events == 1]);
        
        % Do a wilcoxon rank sum test
        p = ranksum(nr_noevent,nr_event,'alpha',0.05,'tail','left');

        ranksum_p(rank_counter) = p;
        rank_counter = rank_counter + 1;
        
    end
    
end

% Find the % of cells that have p < 0.05 and % of cells with p value < 0.01
p_005 = (sum(ranksum_p<0.05)/length(ranksum_p))*100;
p_001 = (sum(ranksum_p<0.01)/length(ranksum_p))*100;


%% Plot
figure(1);clf;

% Show example of pupil image
seqfile = '/Users/annachristinagarvert/Dropbox (UIO Physiology Dropbox)/Lab Data/Andreas Lande/WHISKER CAM DATA/PupilVideo/m1445-20211019-01_pupileTracking_10-15-00.000.seq';
pupil_image = sb.database.load.getSelectedFrames(seqfile,1000);
pupil_image = pupil_image(1:200,1:200);
imagesc(pupil_image);
axis('equal');
axis('off');
colormap('gray');

% Show median and 25th/75th percentile for X and Y and Radius where all mice are pooled
% Create rastermaps of the signals
for s = 1:3
    if s == 1
        rmap_radius = sb.generate.rasterMaps(mData(s).sData.behavior.wheelPosDsBinned,mData(s).sData.behavior.wheelLapDsBinned,mData(s).sData.pupil.radius');
        rmap_x = sb.generate.rasterMaps(mData(s).sData.behavior.wheelPosDsBinned,mData(s).sData.behavior.wheelLapDsBinned,mData(s).sData.pupil.x');
        rmap_y = sb.generate.rasterMaps(mData(s).sData.behavior.wheelPosDsBinned,mData(s).sData.behavior.wheelLapDsBinned,mData(s).sData.pupil.y');
    else
        rmap_radius = cat(1,rmap_radius,sb.generate.rasterMaps(mData(s).sData.behavior.wheelPosDsBinned,mData(s).sData.behavior.wheelLapDsBinned,mData(s).sData.pupil.radius'));
        rmap_x = cat(1,rmap_x,sb.generate.rasterMaps(mData(s).sData.behavior.wheelPosDsBinned,mData(s).sData.behavior.wheelLapDsBinned,mData(s).sData.pupil.x'));
        rmap_y = cat(1,rmap_y,sb.generate.rasterMaps(mData(s).sData.behavior.wheelPosDsBinned,mData(s).sData.behavior.wheelLapDsBinned,mData(s).sData.pupil.y'));
    end
end

% Normalize each lap to itself
for x = 1:size(rmap_radius,1)
   rmap_radius(x,:) = rmap_radius(x,:) - nanmean(rmap_radius(x,:));
   rmap_x(x,:) = rmap_x(x,:) - nanmean(rmap_x(x,:));
   rmap_y(x,:) = rmap_y(x,:) - nanmean(rmap_y(x,:));
end

figure(2); clf;

% Plot the mediam and percentiles of pupil signals along the track
x = 1:105;

% Horizontal
subplot(2,2,1);
plotQuartiles(x,rmap_x(1:50,:),'k',1);
title('Horizontal');
xticks([0,105])
xticklabels({0,157})

% Vertical
subplot(2,2,2);
plotQuartiles(x,rmap_y(1:50,:),'k',1);
title('Vertical');
xticks([0,105])
xticklabels({0,157})

% Radius
subplot(2,2,3);
plotQuartiles(x,rmap_radius(1:50,:),'k',1);
title('Radius');
xticks([0,105])
xticklabels({0,157})

% Add lines to each subplot
for p = 1:3
    subplot(2,2,p);
    hold on;
    ylim1 = get(gca,'YLim');
    line([32,32],[-1,1],'Color','k')
    line([72,72],[-1,1],'Color','k')
    set(gca,'YLim',ylim1,'FontSize',16);
    set(gca,'XLim',[0,105])
end

% P-values of test on radius
subplot(2,2,4);
bar(1:2,[p_005,p_001])
xticklabels({'p < 0.05', 'p < 0.01'})
ylabel('% of landmark cells')
title('Cells correlated with eye radius')
ylim([0,5])
set(gca,'FontSize',16);

set(gcf,'Renderer', 'painters', 'Position', [2000 400 700 600])


