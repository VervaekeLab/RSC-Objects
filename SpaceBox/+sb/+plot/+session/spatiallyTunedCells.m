function spatiallyTunedCells(sData,varargin)
%% spatiallyTunedCells
% Show all place cells found in the session.
%
% INPUT
%   -- Required
%       sData.
%
%   -- Optional
%   varargin: Alternating argument name and value, ex: ["keyword", value].
%       Arguments - See "Default parameters" section.
%
% Written by Andreas S Lande 2019


%% Default parameters
params = struct();
params.spatially_cells_included = {'pcs','tcs','lcs'};

% Sort all rastermaps based on the motor manipulations
try
    params.manipulated_motor_position = sData.landmarks.motor.manipulated; % This can be overwritten as input
catch
    params.manipulated_motor_position = 0;
end

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

mData = struct();
mData(1).sData = sData;
data.n_lcs = [];
data.n_tcs = [];
data.n_pcs = [];
data.n_cells = [];
data.n_acs = [];
data.n_ocs = [];


for s = 1
    
    % Clean up and grab sData for current session
    clearvars -except mData data s params
    sData = mData(s).sData;
    
    % Get data from sData
    pos = sData.behavior.wheelPosDsBinned;
    lap = sData.behavior.wheelLapDsBinned;
    dff = sData.imdata.roiSignals(2).dffSubtractedNpil;
    deconv = smoothdata(sData.imdata.roiSignals(2).deconv','gaussian',10)'; % Smooth the deconvolved signal for a 300 ms (-ish) window
    run = sData.behavior.runSpeedDs;
    
    % Create struct to store all rastermaps
    rmaps = struct();
    
    % Create a binary array of which sample is above or below a running speed threshold
    run_thresholded = [run<2];
    
    % Createx rastermap for deconvolved signal
    rmaps.deconv = sb.generate.rasterMaps(pos,lap,deconv,...
        'smoothing',5, ...
        'exclude_inds', run_thresholded);
    
    % Create rastermap for dF/F signals
    rmaps.dff = sb.generate.rasterMaps(pos,lap,dff, 'exclude_inds', run_thresholded);
    
    % Create rastermap for run speed signal
    rmaps.run = sb.generate.rasterMaps(pos,lap,run,'smoothing',0,'exclude_inds', run_thresholded);
    
    % Create rastermap using the speed signal from the landmark motors
    try
        right_motor = sData.daqdata.right_motor_speed(sData.daqdata.frameIndex);
        rmaps.motor_right = sb.generate.rasterMaps(pos,lap,right_motor,'smoothing',0);
        
        left_motor = sData.daqdata.left_motor_speed(sData.daqdata.frameIndex);
        rmaps.motor_left = sb.generate.rasterMaps(pos,lap,left_motor,'smoothing',0);
        
    catch
        sData = pipe.loadWhiskerTouches(sData);
        sData = sb.helper.convertWhiskerImagingToLandmarkOnsetVector(sData);
    end
    
    % Sort all rastermaps based on the motor manipulations
    manipulated_motor_position = 1; % This means the first motor position is manipulated.
    
    % Divide rmaps into motor on and motor off laps
    laps_on = find([sData.landmarks.motor.on(:,manipulated_motor_position) == 1]);
    laps_off = find([sData.landmarks.motor.on(:,manipulated_motor_position) == 0]);
    
    rmaps.deconv_motor_on = rmaps.deconv(laps_on,:,:);
    rmaps.deconv_motor_off = rmaps.deconv(laps_off,:,:);
    rmaps.dff_motor_on = rmaps.dff(laps_on,:,:);
    rmaps.dff_motor_off = rmaps.dff(laps_off,:,:);
    rmaps.run_motor_on = rmaps.run(laps_on,:);
    rmaps.run_motor_off = rmaps.run(laps_off,:);
    
    % Number of rois
    n_rois = size(rmaps.deconv_motor_on,3);
    
    % Find active cells
    acs = sb.classify.activeRoisAPT(dff,'activity_peak_threshold',0.3);
    
    % Run landmark cells detection algo
    if sum(contains(params.spatially_cells_included,'tcs'))
        [lcs,tcs] = sb.classify.landmarkAndTraceCells(rmaps.deconv_motor_on,rmaps.deconv_motor_off);
        lcs = lcs(ismember(lcs,tcs)==0); % Exclude trace cells from the landmark cells list
    else
        [lcs] = sb.classify.landmarkAndTraceCells(rmaps.deconv_motor_on,rmaps.deconv_motor_on);
        tcs = [];
    end
    
    % Include only cells that show dFF activity
    lcs = lcs(ismember(lcs,find(acs))==1);
    tcs = tcs(ismember(tcs,find(acs))==1);
    
    % Find place cells
    pcs = sb.classify.placeTunedRoisSimple(rmaps.deconv_motor_on);
    pcs = [pcs & acs]; % Include only cells that show dFF activity
    pcs = find(pcs);
    
    % Remove landmark cells from the place cell list so they are unique
    pcs = pcs(ismember(pcs,lcs)==0);
    pcs = pcs(ismember(pcs,tcs)==0);
    
    % Other cells not tuned
    ocs = ones(1,n_rois);
    ocs(unique([lcs,pcs,tcs])) = 0;
    ocs = find(ocs);
    
    % Add data to data
    data.n_lcs = [data.n_lcs,length(lcs)];
    data.n_tcs = [data.n_tcs,length(tcs)];
    data.n_pcs = [data.n_pcs,length(pcs)];
    data.n_cells = [data.n_cells,n_rois];
    data.n_acs = [data.n_acs, sum(acs)];
    data.n_ocs = [data.n_ocs,length(ocs)];
    
    % Add data to mData for later plotting
    mData(s).tcs.deconv_motor_on = rmaps.deconv_motor_on(:,:,tcs);
    mData(s).tcs.deconv_motor_off = rmaps.deconv_motor_off(:,:,tcs);
    mData(s).tcs.dff_motor_on = rmaps.dff_motor_on(:,:,tcs);
    mData(s).tcs.dff_motor_off = rmaps.dff_motor_off(:,:,tcs);
    mData(s).tcs.ind = tcs;
    
    mData(s).lcs.deconv_motor_on = rmaps.deconv_motor_on(:,:,lcs);
    mData(s).lcs.deconv_motor_off = rmaps.deconv_motor_off(:,:,lcs);
    mData(s).lcs.dff_motor_on = rmaps.dff_motor_on(:,:,lcs);
    mData(s).lcs.dff_motor_off = rmaps.dff_motor_off(:,:,lcs);
    mData(s).lcs.ind = lcs;
    
    lcs_not_tcs = lcs(ismember(lcs,tcs)==0);
    mData(s).lcs_not_tcs.deconv_motor_on = rmaps.deconv_motor_on(:,:,lcs_not_tcs);
    mData(s).lcs_not_tcs.deconv_motor_off = rmaps.deconv_motor_off(:,:,lcs_not_tcs);
    mData(s).lcs_not_tcs.ind = lcs_not_tcs;
    mData(s).pcs.deconv_motor_on = rmaps.deconv_motor_on(:,:,pcs);
    mData(s).pcs.deconv_motor_off = rmaps.deconv_motor_off(:,:,pcs);
    mData(s).pcs.dff_motor_on = rmaps.dff_motor_on(:,:,pcs);
    mData(s).pcs.dff_motor_off = rmaps.dff_motor_off(:,:,pcs);
    mData(s).pcs.ind = pcs;
    
    mData(s).ocs.deconv_motor_on = rmaps.deconv_motor_on(:,:,ocs);
    mData(s).ocs.deconv_motor_off = rmaps.deconv_motor_off(:,:,ocs);
    mData(s).ocs.dff_motor_on = rmaps.dff_motor_on(:,:,ocs);
    mData(s).ocs.dff_motor_off = rmaps.dff_motor_off(:,:,ocs);
    mData(s).ocs.ind = ocs;
    
    
    % Visualize all cells
    %sb.plot.rasterMaps({rmaps.deconv_motor_on(:,:,tcs),rmaps.deconv_motor_off(:,:,tcs)},'rmap_title',num2cell(tcs))
    
end

% Create raster map titles
s = 1;

% Plot place cells
if sum(contains(params.spatially_cells_included,'pcs'))>0
    rmaptitles = {};
    for i = 1:length(mData(s).pcs.ind)
        rmaptitles{i} = sprintf('pc %i',mData(s).pcs.ind(i));
    end

    sb.plot.rasterMaps({mData(s).pcs.dff_motor_on,mData(s).pcs.dff_motor_off},'rmap_title',rmaptitles','add_lines_vertical',[33,73,93]);
end

% Plot landmark cells
if sum(contains(params.spatially_cells_included,'lcs'))>0
    rmaptitles = {};
    for i = 1:length(mData(s).lcs.ind)
        rmaptitles{i} = sprintf('lc %i',mData(s).lcs.ind(i));
    end

    sb.plot.rasterMaps({mData(s).lcs.dff_motor_on,mData(s).lcs.dff_motor_off},'rmap_title',rmaptitles,'add_lines_vertical',[33,73,93]);
end

% Plot trace cells
if sum(contains(params.spatially_cells_included,'tcs'))>0
    rmaptitles = {};
    for i = 1:length(mData(s).tcs.ind)
        rmaptitles{i} = sprintf('tc %i',mData(s).tcs.ind(i));
    end

    sb.plot.rasterMaps({mData(s).tcs.dff_motor_on,mData(s).tcs.dff_motor_off},'rmap_title',rmaptitles,'add_lines_vertical',[33,73,93]);
end


end