% parameters
window_size = 31; % time bins before and after the landmark interaction
z_score_threshold = 4;
landmark_color = [239,35,60]/255;

%% FIGURE B:  Plot 3 example predictive cells
figure(2); clf;
example_cells = [20,14,50];
axis_grid = reshape(1:21,[3,7])';

% For each of the example cells
for c = 1:3
    
    % Get correct session number and roi
    s = predictive_cells_session(example_cells(c))
    roi = predictive_cells_inds(example_cells(c));
    
    % Create a rastermap of the deconv responses around the landmark at the
    % first landmark position, where x is time and y is lap
    deconv_at_landmark = getDeconvSignalAroundLandmark(mData(s).sData,roi);
    
    % Draw rastermap of the predictive cell
    subplot(7,3,axis_grid(1:end-4,c));
    rmap = predictive_rmaps(1:50,:,example_cells(c));
    rmap = normalize(rmap);
%     rmap = zScoreNormalize(rmap,'all');
    imagesc(rmap,[0,z_score_threshold]);
    hold on
    
    % Set lines for indicating where the landmark motor starts to move.
    % I.e. there is no touch of the whiskers at this point.
    line([30,30],[0,100],'Color',landmark_color,'LineWidth',0.5);
    line([70,70],[0,100],'Color',landmark_color,'LineWidth',0.5);
    
    % Plot the average response at the bottom
    plot(50-normalize(nanmean(rmap),'range',[0,10]),'LineWidth',1,'Color','k');
        
    % Properties
    if c == 1
        ylabel('Lap #')
    end
    xticks([]);
    yticks([1,50])
    title(sprintf('Cell #%i',c))
    set(gca,'FontSize',18);
    colormap(flipud(colormap('gray')));
    
    % Draw average
    subplot(7,3,axis_grid(end-3:end-2,c));
    
    % Create a rastermap of the running speed and use this to compute
    % percentiles for a shaded error plot
    run_rmap = createRunningSpeedRasterMap(mData(s).sData);
    
    % Compute average running speed and the max speed
    average_run = nanmean(run_rmap(1:50,:));
    max_run_speed = round(max(average_run));
    
    
    
    plot(nanmedian(run_rmap),'Color','k');
    hold on;
    
    % Plot 25th percentile and 75th percentile
    plot(prctile(run_rmap,25),'Color','b');
    plot(prctile(run_rmap,75),'Color','b');
    
    % Set lines for indicating where the landmark motor starts to move.
    % I.e. there is no touch of the whiskers at this point.
    line([30,30],[-2,200],'LineWidth',0.5,'Color',landmark_color)
    line([70,70],[-2,200],'LineWidth',0.5,'Color',landmark_color)
    xlim([0,105]);
    yticks([0,max_run_speed]);
    
    if c == 1
        ylabel('Speed (cm/s)')
    end
    
    ylim([0,max_run_speed+max_run_speed*0.2])
    xticks([1,105])
    xticklabels({1,157})
    xlabel('Position (cm)')
    set(gca,'FontSize',16)
        
    
    % Draw deconvolved response before landmark onset (time axis)
    subplot(7,3,axis_grid(end-1:end,c));
    plot(normalize(nanmean(smoothdata(deconv_at_landmark,2,'movmean',5)),'range',[0,1]),'LineWidth',1,'Color','k');
    hold on;
    line([window_size,window_size],[-2,2],'LineWidth',0.5,'Color',landmark_color)
    yticks([0,1]);
    ylim([-0.05,1.05])
    xticks([1,window_size,window_size*1.5, window_size*2]);
    xticklabels({-1,0,0.5})
    xlim([0,window_size*1.5])
    
    if c == 1
       ylabel('Norm. deconv'); 
    end
    
    xlabel('Time (sec)')
    set(gca,'FontSize',16)
    
    
    
end

% Set correct size of the figure
set(gcf,'renderer', 'painters', 'Position', [2200,100,500,600])



%% --- HELPER FUNCTIONS

function rmap = getDeconvSignalAroundLandmark(sData,roi)
% Returns a rastermap where each row is the deconvolved signal n bins
% before and after the position where the motor begins to move at the first
% landmark position.

% Parameters
landmark_position = 48; % This is the position where the motor begins to move.
window_size = 31; % time bins before and after the landmark interaction

% Get data to use
deconv = sData.imdata.roiSignals(2).deconv(roi,:);
pos = sData.behavior.wheelPosDs;

% Set some init vars
bins_at_landmark = [];
landmark_passed = true;
i = 0;

% Find each bin where the mouse enters position 45 cm.
while i < length(pos)
    i = i+1;
    if pos(i) > landmark_position

        if ~landmark_passed
            bins_at_landmark = [bins_at_landmark,i];
            landmark_passed = true;
            i = i+20; % Do a jump so that we don't get multiple values in case there are more interactions with the landmark
        end
    else
        landmark_passed = false;
    end

end

% Get the deconv signal at each landmark passing
window_n_bins = length([1-window_size:1+window_size]);
n_laps_to_include = 50;

rmap = nan(n_laps_to_include,window_n_bins);
for i = 1:n_laps_to_include
    rmap(i,:) = deconv(bins_at_landmark(i)-window_size:bins_at_landmark(i)+window_size);
end

end


function rmap = createRunningSpeedRasterMap(sData)

% Get data
pos = sData.behavior.wheelPosDsBinned;
lap = sData.behavior.wheelLapDsBinned;
speed = sData.behavior.runSpeedDs;

% Create rmap
rmap = sb.generate.rasterMaps(pos,lap,speed,'smoothing',3);

end




