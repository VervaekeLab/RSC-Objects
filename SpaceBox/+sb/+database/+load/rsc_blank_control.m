function [mData,data] = rsc_blank_control()
%% rsc_blank_control
% Load the recordings where the mouse has had whiskers intact and recorded
% in darkness, but there were no landmarks (nor sandpaper cue) on wheel.

% Recordings to include
data.sessionIDs = { 'm1442_20211004_01', ... 
                    'm1444_20211004_01', ... 
                    'm1444_20211004_02', ... 
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
    mData(s).rmaps = sb.classify.cells(mData(s).sData,'manipulated_motor',0,'dff_signal','dffSubtractedNpil');
    
    % Add cell statistics to data struct
    data.n_cells(end+1) = mData(s).rmaps.meta.n_cells;
    data.n_lcs(end+1) = mData(s).rmaps.meta.n_lcs;
    data.n_tcs(end+1) = mData(s).rmaps.meta.n_tcs;
    data.n_pcs(end+1) = mData(s).rmaps.meta.n_pcs;
    data.n_acs(end+1) = mData(s).rmaps.meta.n_acs;
    data.n_ocs(end+1) = mData(s).rmaps.meta.n_ocs;
    
    % Get mouse name from each recording
    data.mouseNames{s} = data.sessionIDs{s}(1:5);
    
end

% Print text info
fprintf('\nTotal cells: %i \nActive cells: %i\nAll spatially tuned cells: %i (%.1f%%)\nPlace cells: %i (%.1f%%)\nLandmark cells: %i (%.1f%%)\nTrace cells: %i (%.1f%%)\n',sum(data.n_cells),sum(data.n_acs),sum([data.n_pcs,data.n_lcs,data.n_tcs]),(sum([data.n_pcs,data.n_lcs,data.n_tcs])/sum(data.n_acs))*100,sum(data.n_pcs),(sum(data.n_pcs)/sum(data.n_acs))*100,sum(data.n_lcs),(sum(data.n_lcs)/sum(data.n_acs))*100,sum(data.n_tcs),(sum(data.n_tcs)/sum(data.n_acs))*100);
