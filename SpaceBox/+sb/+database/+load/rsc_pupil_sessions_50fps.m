function [mData,data] = rsc_pupil_sessions_50fps()
%% rsc_pupil_sessions
% Load the recordings where the mouse has the whiskers intact. The first 50
% laps of these recordings are in darkness, then next 50 are in 0.03 V
% light, then the remainder is 0.1V light. 

% These recordings are all at 50 fps.

% Recordings to include
data.sessionIDs = { 'm1445_20211019_01', ...
                    'm1446_20211019_01', ...
                    'm1447_20211019_01', ...
                    };

% Load all sessions
mData = sb.database.loadSessions(data.sessionIDs);

%% Compute rastermaps and classify cells for each session

% Initiate some variables
data.n_cells = [];
data.n_lcs = [];
data.n_tcs = [];
data.n_pcs = [];
data.n_acs = [];
data.n_ocs = [];

% Create rastermaps for each session
for s = 1:length(mData)
    
% Downsample the light voltage signal
    mData(s).sData.daqdata.lightDS = mData(s).sData.daqdata.light_voltage(mData(s).sData.daqdata.frameIndex);
    
    % Get data from sData
    sData = mData(s).sData;
    
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
    
    mData(s).rmaps = rmaps;
    
    % Get light voltage signal
    light = mData(s).sData.daqdata.lightDS;
    
    % Create rastermap for the light
    mData(s).rmaps.light.all = sb.generate.rasterMaps(pos,lap,light);
   
    % For each lap, grab average voltage
    light_levels_per_lap = round(nanmean(mData(s).rmaps.light.all,2),3);
    
    laps_light_level_1 = find(light_levels_per_lap == 0);
    laps_light_level_2 = find(light_levels_per_lap == 0.03);
    laps_light_level_3 = find(light_levels_per_lap == 0.1);
    
    mData(s).rmaps.level1.dff = rmaps.dff(laps_light_level_1,:,:);
    mData(s).rmaps.level1.deconv = rmaps.deconv(laps_light_level_1,:,:);
    mData(s).rmaps.level2.dff = rmaps.dff(laps_light_level_2,:,:);
    mData(s).rmaps.level2.deconv = rmaps.deconv(laps_light_level_2,:,:);
    mData(s).rmaps.level3.dff = rmaps.dff(laps_light_level_3,:,:);
    mData(s).rmaps.level3.deconv = rmaps.deconv(laps_light_level_3,:,:);
    
    % Create rastermap of the light signal for control
    mData(s).rmaps.light.level1 = mData(s).rmaps.light.all(laps_light_level_1,:);
    mData(s).rmaps.light.level2 = mData(s).rmaps.light.all(laps_light_level_2,:);
    mData(s).rmaps.light.level3 = mData(s).rmaps.light.all(laps_light_level_3,:);
    
    % Load pupil data
    mData(s).sData = sb.database.loadPupilData(mData(s).sData); % this function is in pipe.  and not sb in andreas set
    
    % Convert pupil data to match the number of bins from calcium imaging.
    imaging_n_frames = size(mData(s).sData.imdata.roiSignals(2).dff,2);
    pupil_n_frames = length(mData(s).sData.pupil.x);
    mData(s).sData.pupil.x = sb.database.load.resampleX(mData(s).sData.pupil.x,imaging_n_frames,pupil_n_frames);
    mData(s).sData.pupil.y = sb.database.load.resampleX(mData(s).sData.pupil.y,imaging_n_frames,pupil_n_frames);
    mData(s).sData.pupil.height = sb.database.load.resampleX(mData(s).sData.pupil.height,imaging_n_frames,pupil_n_frames);
    mData(s).sData.pupil.width = sb.database.load.resampleX(mData(s).sData.pupil.width,imaging_n_frames,pupil_n_frames);
    mData(s).sData.pupil.orientation = sb.database.load.resampleX(mData(s).sData.pupil.orientation,imaging_n_frames,pupil_n_frames);
    mData(s).sData.pupil.radius = sb.database.load.resampleX(mData(s).sData.pupil.radius,imaging_n_frames,pupil_n_frames);
    
end

% Find the number of place cells and landmark cells at each light level
for s = 1:length(mData)
    for level = {'level1','level2','level3'}
        
        level = level{1};
%         mData(s).rmaps = sb.classify.cells(mData(s).sData,data.sessionIDs{s},'manipulated_motor',1,'dff_signal','dffSubtractedNpil');
        
        pcs = sb.classify.placeCellsShuffle(mData(s).rmaps.(level).deconv);
        lcs = sb.classify.landmarkCellsByShuffling(mData(s).rmaps.(level).deconv);
        
        % Remove lcs from pcs list
        pcs = find(pcs([~ismember(find(pcs),lcs)]));
        
        mData(s).rmaps.(level).pcs = pcs;
        mData(s).rmaps.(level).lcs = lcs;
        
    end
    
end

% Print text info
fprintf('\nTotal cells: %i \nActive cells: %i\nAll spatially tuned cells: %i (%.1f%%)\nPlace cells: %i (%.1f%%)\nLandmark cells: %i (%.1f%%)\nTrace cells: %i (%.1f%%)\n',sum(data.n_cells),sum(data.n_acs),sum([data.n_pcs,data.n_lcs,data.n_tcs]),(sum([data.n_pcs,data.n_lcs,data.n_tcs])/sum(data.n_acs))*100,sum(data.n_pcs),(sum(data.n_pcs)/sum(data.n_acs))*100,sum(data.n_lcs),(sum(data.n_lcs)/sum(data.n_acs))*100,sum(data.n_tcs),(sum(data.n_tcs)/sum(data.n_acs))*100);
