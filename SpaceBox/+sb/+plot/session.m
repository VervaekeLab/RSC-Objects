function session(sData,varargin)
%% session
% The function does basic visualizations of a sData session to see that it
% all looks ok.
%
% INPUT
%   -- Required
%       sData: sData struct.
%
% Written by Andreas S Lande 2019

%% Default parameters
params = struct();
params.signal = 'dffSubtractedNpil'; 
params.show_roistats = false;
params.plot_summary = true; % Summary plot of motors, runspeed, dFF / deconv signal etc
params.plot_all_cells = false; % Plots the rastermap for each cell 
params.plot_place_cells = false;
params.plot_trace_cells = false;
params.plot_landmark_cells = false;

% Print all variables if "help" is only input 
try
    if inp == "help"
        global verbose; verbose = true;
        sb.helper.displayParameters(params)
        verbose = false;
        return;
    end
end

% Update parameters
params = sb.helper.updateParameterStruct(params,varargin);


%% Create data to be plotted
if params.show_roistats
    roiStats = getRoiActivityStatsFromSignalArray(sData.imdata.roiSignals(2).(params.signal));
end

pos = sData.behavior.wheelPosDsBinned;
lap = sData.behavior.wheelLapDsBinned;
dff = sData.imdata.roiSignals(2).(params.signal);
deconv = sData.imdata.roiSignals(2).deconv;
run = sData.behavior.runSpeedDs;

% Create dff rastermap
rmaps.dff = sb.generate.rasterMaps(pos,lap,dff);

% Create deconv rastermap
rmaps.deconv = sb.generate.rasterMaps(pos,lap,deconv,'smoothing',10);

% Create motor rastermap
right_motor = sData.daqdata.right_motor_speed(sData.daqdata.frameIndex);
rmaps.motor_right = sb.generate.rasterMaps(pos,lap,right_motor,'smoothing',0);

left_motor = sData.daqdata.left_motor_speed(sData.daqdata.frameIndex);
rmaps.motor_left = sb.generate.rasterMaps(pos,lap,left_motor,'smoothing',0);

% Create run speed map
rmaps.run = sb.generate.rasterMaps(pos,lap,run,'smoothing',0);

% Sort all rastermaps based on the motor manipulations
manipulated_motor_position = 1; % sData.landmarks.motor.manipulated;

% If landmark is manipulated, sort rastermaps thereafter
if manipulated_motor_position
    rmaps.run = sb.sort.rasterMapByLandmarkOnset(rmaps.run,sData.landmarks.motor.on(:,manipulated_motor_position));
    rmaps.motor_right_sorted = sb.sort.rasterMapByLandmarkOnset(rmaps.motor_right,sData.landmarks.motor.on(:,manipulated_motor_position));
    rmaps.motor_left_sorted = sb.sort.rasterMapByLandmarkOnset(rmaps.motor_left,sData.landmarks.motor.on(:,manipulated_motor_position));
end


%% Figure 1
if params.plot_summary
    figure(1);clf;

    subplot(3,5,1);
    imagesc(sData.imdata.meta.fovImageMax);
    colormap('gray');
    title('FOV image');
    xticks([]);
    yticks([]);
    set(gca,'FontSize',15)

    if params.show_roistats
        subplot(3,5,2);
        violin(roiStats.peakDff,'facecolor',[1,0,0])
        title('Peak dF/F');
        xticks([]);
        set(gca,'FontSize',15)

        subplot(3,5,6);
        violin(roiStats.signalToNoise,'facecolor',[0,1,0])
        title('Signal to noise');
        xticks([]);
        set(gca,'FontSize',15)
        ylabel(sprintf('SessionID: %s',sData.sessionInfo.sessionID(1:17)),'FontSize',25);

        subplot(3,5,7);
        violin(roiStats.activityLevel,'facecolor',[0,0,1])
        title('Activity level');
        ylim([0,1]);
        xticks([]);
        set(gca,'FontSize',15)
    end

    subplot(3,5,[5,10]);
    imagesc(rmaps.run);
    title('Running speed');
    set(gca,'FontSize',15)
    ylabel('Lap #');
    % Create corrext xtix based on bin size
    xtix = get(gca, 'XTick');
    set(gca, 'XTick',xtix, 'XTickLabel',xtix*1.5);

    subplot(3,5,4);
    imagesc(rmaps.motor_right);
    title('Right Motor');
    set(gca,'FontSize',15)
    % Create corrext xtix based on bin size
    xtix = get(gca, 'XTick');
    set(gca, 'XTick',xtix, 'XTickLabel',xtix*1.5);

    subplot(3,5,3);
    imagesc(rmaps.motor_left);
    title('Left Motor');
    set(gca,'FontSize',15)
    % Create corrext xtix based on bin size
    xtix = get(gca, 'XTick');
    set(gca, 'XTick',xtix, 'XTickLabel',xtix*1.5);

    subplot(3,5,8);
    imagesc(rmaps.motor_left_sorted);
    title('Left Motor Sorted');
    set(gca,'FontSize',15)
    % Create corrext xtix based on bin size
    xtix = get(gca, 'XTick');
    set(gca, 'XTick',xtix, 'XTickLabel',xtix*1.5);

    subplot(3,5,9);
    imagesc(rmaps.motor_right_sorted);
    title('Right Motor Sorted');
    set(gca,'FontSize',15)
    % Create corrext xtix based on bin size
    xtix = get(gca, 'XTick');
    set(gca, 'XTick',xtix, 'XTickLabel',xtix*1.5);

    % Find all place cells and sort their peaks
    pcs = sb.classify.placeTunedRoisSimple(rmaps.deconv);
    pcs = find(pcs);
    rmaps.sorted_place_cells_deconv = sb.sort.rasterMapByPeak(sb.generate.averagedTuningFromRasterMaps(rmaps.deconv(:,:,pcs)));
    rmaps.sorted_place_cells_dff = sb.sort.rasterMapByPeak(sb.generate.averagedTuningFromRasterMaps(rmaps.dff(:,:,pcs)));


    subplot(3,5,14);
    imagesc(rmaps.sorted_place_cells_dff);
    title(sprintf('dF/F - %.1f %% (%i/%i)',(length(pcs)/size(rmaps.deconv,3))*100,length(pcs),size(rmaps.deconv,3)));
    set(gca,'FontSize',15)
    % Create corrext xtix based on bin size
    xtix = get(gca, 'XTick');
    set(gca, 'XTick',xtix, 'XTickLabel',xtix*1.5);
    ylabel('Cell #');
    xlabel('Position (cm)');

    subplot(3,5,15);
    imagesc(rmaps.sorted_place_cells_deconv);
    title(sprintf('deconv - %.1f %% (%i/%i)',(length(pcs)/size(rmaps.deconv,3))*100,length(pcs),size(rmaps.deconv,3)));
    set(gca,'FontSize',15)
    % Create corrext xtix based on bin size
    xtix = get(gca, 'XTick');
    set(gca, 'XTick',xtix, 'XTickLabel',xtix*1.5);
    ylabel('Cell #');
    xlabel('Position (cm)');

    % Plot example place cell activity and position on wheel
    subplot(3,5,11:13);
    cla;
    end_ind = 3000;
    %plot(normalize(sData.behavior.wheelPosDsBinned(1:end_ind),'range'),'LineWidth',3,'Color','k');
    ylim([-0.1,4.1]);
    hold on
    shift = 0.1;
    for x = 1:10

        plot(normalize(sData.imdata.roiSignals(2).dff(pcs(x),1:end_ind),'range')+shift);
        hold on
        plot(normalize(sData.imdata.roiSignals(2).deconv(pcs(x),1:end_ind),'range')+shift,'Color','k');
        shift = shift + 0.4;

    end
end

% Figure 2
if params.plot_all_cells
    % TODO
   rmaps = rmapsCreate(sData,'all');
    
   sb.plot.rasterMaps(rmaps);
    
end

% Figure 3
if params.plot_place_cells
    % TODO
    % Plot all place cells
    error('Not implemented')
end

% Figure 4
if params.plot_landmark_cells
    % TODO
    error('Not implemented')
end

% Figure 5
if params.plot_trace_cells
    % TODO
    error('Not implemented')
end

end

% ---- HELPER FUNCTIONS
function rmaps = rmapsCreate(sData,ctype)

    pos = sData.behavior.wheelPosDsBinned;
    laps = sData.behavior.wheelLapDsBinned;
    run = sData.behavior.runSpeedDs;
    signal = sData.imdata.roiSignals(2).dffSubtractedNpil;
    
    rmaps = sb.generate.rasterMaps(pos,laps,signal,'exclude_inds',run<2);
    
end