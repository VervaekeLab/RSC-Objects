
function [mData,data] = rsc_introduction(data,mData)
  %% analysis of tactile object introduction
%
% INPUT
%   -- Required: dirct
%
% OUTPUT
%      mData containing rmaps
%      data with classification
%

% set seed to make reproducible
% please add loading direction to sData structs here!
dirct = '';

data.sessionIDs = { 'm1551-20230227-03', ... # Landmarks removed after 41 laps. Ran 96 laps total
           'm1551-20230529-01',...
            'm1553-20230227-01', ... # Landmarks removed after 41 laps. Ran 136 laps total
            'm1554-20230227-01', ... # Landmarks removed after 30 laps. Ran 44 laps total
            'm1554-20230529-01',...
            'm7131-20230529-01',...
            'm7132-20230529-01'};
    

% Load all sessions

for f = 1:length(data.sessionIDs)
    
    load(fullfile(dirct,data.sessionIDs{f},'sData.mat'));
    mData(f).sData = sData;
end



% Initiate some variables
data.n_cells_on = [];
data.n_lcs_on = [];
data.n_tcs_on = [];
data.n_pcs_on = [];
data.n_acs_on = [];
data.n_ocs_on = [];

data.n_cells_off = [];
data.n_lcs_off = [];
data.n_tcs_off = [];
data.n_pcs_off = [];
data.n_acs_off = [];
data.n_ocs_off = [];

% Create rastermaps for each session
for s = 1:length(mData)
    
    mData(s).sData.imdata.roiSignals(2).('dffSubtractedNpil') = mData(s).sData.imdata.roiSignals(2).dff;
    rmaps = sb.generate.rasterMapsForSession(mData(s).sData,"dff_signal", "dffSubtractedNpil",'manipulated_motor',1);

    on_laps = find([mData(s).sData.landmarks.motor.on(:,1) == 1]);
    off_laps = find([mData(s).sData.landmarks.motor.on(:,1) == 0]);

    % identify landmark and place cells in on trials 
    mData(s).rmaps_on = sb.classify.classifyCells(mData(s).sData,rmaps.deconv(on_laps,:,:));  
    % Add cell statistics to data struct 
    data.n_cells_on(end+1) = mData(s).rmaps_on.meta.n_cells;
    data.n_lcs_on(end+1)   = mData(s).rmaps_on.meta.n_lcs;
    data.n_pcs_on(end+1)   = mData(s).rmaps_on.meta.n_pcs;
    data.n_acs_on(end+1)   = size(rmaps.deconv,3);%mData(s).rmaps_on.meta.n_acs;
    data.n_ocs_on(end+1)   = mData(s).rmaps_on.meta.n_ocs;
    mData(s).rmaps_on.pcs.deconv_off = rmaps.deconv(off_laps,:,mData(s).rmaps_on.pcs.ind);
    mData(s).rmaps_on.lcs.deconv_off = rmaps.deconv(off_laps,:,mData(s).rmaps_on.lcs.ind);

     % identify landmark and place cells in on trials 
     mData(s).rmaps_off = sb.classify.classifyCells(mData(s).sData,rmaps.deconv(off_laps,:,:));  
    % Add cell statistics to data struct 
    data.n_cells_off(end+1) = mData(s).rmaps_off.meta.n_cells;
    data.n_lcs_off(end+1)   = mData(s).rmaps_off.meta.n_lcs;
    data.n_pcs_off(end+1)   = mData(s).rmaps_off.meta.n_pcs;
    data.n_acs_off(end+1)   = size(rmaps.deconv,3);%mData(s).rmaps_off.meta.n_acs;
    data.n_ocs_off(end+1)   = mData(s).rmaps_off.meta.n_ocs;
    mData(s).rmaps_off.pcs.deconv_on = rmaps.deconv(off_laps,:,mData(s).rmaps_off.pcs.ind);
    mData(s).rmaps_off.lcs.deconv_on = rmaps.deconv(off_laps,:,mData(s).rmaps_off.lcs.ind);
       
end

end


