%% Plot the dF/F and deconvolved signal for four example landmark cells, add position and running speed.

% Parameters
session = 1;
lineWidth = 0.5;
deconv_color = [0,0,250]/255;%[251,176,59]/255;
ind_start = 2270; % bin in recording
ind_end = 2950; % bin in recording

% Plot the landmarks
figure(1); clf;
c = 1; % init

% For each of the selected landmark cells
for cell_index = [3,17]  %[3,5,4,13]   tcs: [3,9,17], pcs: [11,45]
    
    c = c+1;
    
    sData = mData(session).sData;
    
    % Plot activity
    dffToPlot = normalize(sData.imdata.roiSignals(2).dffSubtractedNpil(mData(session).rmaps.tcs.ind(cell_index),ind_start:ind_end),'range',[0,1]);
    deconvToPlot = smoothdata(-0.2+normalize(sData.imdata.roiSignals(2).deconv(mData(session).rmaps.tcs.ind(cell_index),ind_start:ind_end),'range',[0,0.5]),'gaussian',3);

    deconvToPlot = deconvToPlot + c;
    dffToPlot = dffToPlot + c;
    
    plot(dffToPlot,'LineWidth',lineWidth,'Color','k')
    hold on
    plot(deconvToPlot,'LineWidth',lineWidth,'Color',deconv_color)
    xlim([0,ind_end-ind_start])
    ylim([0,6.2])
    yticks([])
    xticks([])
    
    box off
end

% For each of the selected  place cells
for cell_index = [11,45]
    
    c = c+1;
    
    sData = mData(session).sData;
    
    % Plot activity
    dffToPlot = normalize(sData.imdata.roiSignals(2).dffSubtractedNpil(mData(session).rmaps.pcs.ind(cell_index),ind_start:ind_end),'range',[0,1]);
    deconvToPlot = smoothdata(-0.2+normalize(sData.imdata.roiSignals(2).deconv(mData(session).rmaps.pcs.ind(cell_index),ind_start:ind_end),'range',[0,0.5]),'gaussian',3);

    deconvToPlot = deconvToPlot + c;
    dffToPlot = dffToPlot + c;
    
    plot(dffToPlot,'LineWidth',lineWidth,'Color','k')
    hold on
    plot(deconvToPlot,'LineWidth',lineWidth,'Color',deconv_color)
    xlim([0,ind_end-ind_start])
    ylim([0,6.2])
    yticks([])
    xticks([])
    
    box off
end


% Plot position
plot(normalize(sData.behavior.wheelPosDsBinned(ind_start:ind_end),'range',[1,1.8]),'LineWidth',0.75,'Color','k');


% Plot running speed
plot(normalize(sData.behavior.runSpeedDs(ind_start:ind_end),'range',[0.1,0.8]),'LineWidth',0.75,'Color','k');
max_run_speed = round(max(sData.behavior.runSpeedDs(ind_start:ind_end)));

% Add labels to y axis
yticks([0.1,0.8,1,1.8])
yticklabels({0,max_run_speed,0,157})

% Draw a line to indicate 5 seconds
line([100, 255],[6,6],'Color','k','LineWidth',1);

% Set figure size
set(gcf,'Renderer', 'painters', 'Position', [1200 200 700 900])
set(gca,'FontSize',16);

%% Plot FOV
figure;clf;
imagesc(mData(session).sData.imdata.meta.fovImageMax)
colormap('gray')
xticks([]);
yticks([]);
