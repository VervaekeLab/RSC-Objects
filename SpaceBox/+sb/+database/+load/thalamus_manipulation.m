function [mData,data] = thalamus_manipulation()
%% thalamus_manipulation
% Load the dataset for the manipulation experiments in Thalamus. Do basic
% processing.

% Recordings to include
data.sessionIDs = { 'm1807_20210630_01', ...
                    'm1807_20210701_01', ...
                    'm1808_20210701_01', ...
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
fprintf('\nTotal cells: %i \nActive cells: %i\nAll spatially tuned cells: %i (%.1f%%)\nPlace cells: %i (%.1f%%)\nLandmark cells: %i (%.1f%%)\nTrace cells: %i (%.1f%%)\n',sum(data.n_cells),sum(data.n_acs),sum([data.n_pcs,data.n_lcs,data.n_tcs]),(sum([data.n_pcs,data.n_lcs,data.n_tcs])/sum(data.n_acs))*100,sum(data.n_pcs),(sum(data.n_pcs)/sum(data.n_acs))*100,sum(data.n_lcs),(sum(data.n_lcs)/sum(data.n_acs))*100,sum(data.n_tcs),(sum(data.n_tcs)/sum(data.n_acs))*100);
