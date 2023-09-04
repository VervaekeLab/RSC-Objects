function [mData,data] = rsc_manipulation()
%% rsc_manipulation
% Load the dataset for the manipulation experiments in RSC. Do basic
% processing.

% Recordings to include
data.sessionIDs = { 'm1410_20191023_01', ... # 1 OLD (10-14 ish days training with cue setup)
                    'm1411_20191023_01', ... # 2 OLD (10-14 ish days training with cue setup)
                    'm1412_20191023_01', ... # 3 OLD (10-14 ish days training with cue setup)
                    'm1415_20191023_01', ... # 4 OLD (11 ish days training with cue setup)
                    'm1420_20200225_01', ... # 5 NEW (125 ym depth) (15 days training with cue setup)
                    'm1420_20200304_01', ... # 6 NEW (180 ym depth) (20 days training with cue setup)
                    'm1422_20200225_01', ... # 7 NEW (130 ym depth) (9 days training with cue setup)
                    'm1424_20200311_01', ... # 8 NEW (Posterior FOV) (6 days training with cue setup)
                    'm1424_20200315_01', ... # 9 NEW (Anterior FOV) (9 days training with cue setup) %'m1425_20200311_01', ... # 10 NEW (11 days training on wheel)
                    'm1438_20210606_01', ... # 11 NEW (5 days training with landmarks)
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
    
    mData(s).rmaps = sb.classify.cells(mData(s).sData,'manipulated_motor',1,'dff_signal','dffSubtractedNpil');
    
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
fprintf('\nTotal cells: %i \nActive cells: %i\nAll spatially tuned cells: %i (%.1f%%)\nPlace cells: %i (%.1f%%)\nPure landmark cells: %i (%.1f%%)\nTrace cells: %i (%.1f%%)\nAll landmark cells: %i (%.1f%%)\n',sum(data.n_cells),sum(data.n_acs),sum([data.n_pcs,data.n_lcs,data.n_tcs]),(sum([data.n_pcs,data.n_lcs,data.n_tcs])/sum(data.n_acs))*100,sum(data.n_pcs),(sum(data.n_pcs)/sum(data.n_acs))*100,sum(data.n_lcs),(sum(data.n_lcs)/sum(data.n_acs))*100,sum(data.n_tcs),(sum(data.n_tcs)/sum(data.n_acs))*100,sum([data.n_tcs,data.n_lcs]),(sum([data.n_tcs,data.n_lcs])/sum(data.n_acs))*100);
