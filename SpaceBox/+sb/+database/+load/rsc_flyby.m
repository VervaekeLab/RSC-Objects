function [mData,data] = rsc_flyby()
%% rsc_flyby
% Load the dataset for the flyby experiments in RSC. Do basic
% processing.

% Recordings to include
data.sessionIDs = { %'m1420_20200315_01', ...
                    %'m1420_20200315_03', ...
                    'm5115_20210222_02', ...
                    'm5115_20210222_04', ... % this session do not have moving flyby
                    'm5115_20210225_01', ...
                    'm5115_20210225_02' ...
                    'm1438_20210606_01', ...
                    'm1438_20210606_02', ...
                    'm1446_20210930_01', ...
                    'm1446_20210930_02', ...
                    };

% Initiate some variables
data.n_cells = [];
data.n_lcs = [];
data.n_tcs = [];
data.n_pcs = [];
data.n_acs = [];
data.n_ocs = [];

% Load all sessions
mData = sb.database.loadSessions(data.sessionIDs);

%% Compute rastermaps for cells that exist in both of the manipulation session and flyby session pairs
% Filter out cells that are not present in both of paired sessions
for s = 1:2:length(mData)-1
    [mData(s).sData, mData(s+1).sData] = sb.helper.filterOutUniqueRoiIDsAcrossTwoSessions(mData(s).sData,mData(s+1).sData);

    % For each recording find landmark / trace cells
    mData(s).rmaps = sb.classify.cells(mData(s).sData,data.sessionIDs{s},'manipulated_motor',1);

    % Plot the landmark cells and their activity when flyby occurs
    mData(s+1).rmaps.flyby.dff = sb.generate.flyByRasterMaps(mData(s+1).sData);
    mData(s+1).rmaps.flyby.dffZ = sb.generate.flyByRasterMaps(mData(s+1).sData,'signal_type','dffZScored');
    mData(s+1).rmaps.flyby.deconv = sb.generate.flyByRasterMaps(mData(s+1).sData,'signal_type','deconv');
    
    % Add cell statistics to data struct
    data.n_cells(end+1) = mData(s).rmaps.meta.n_cells;
    data.n_lcs(end+1) = mData(s).rmaps.meta.n_lcs;
    data.n_tcs(end+1) = mData(s).rmaps.meta.n_tcs;
    data.n_pcs(end+1) = mData(s).rmaps.meta.n_pcs;
    data.n_acs(end+1) = mData(s).rmaps.meta.n_acs;
    data.n_ocs(end+1) = mData(s).rmaps.meta.n_ocs;
    
end

% Print text info
fprintf('\nTotal cells: %i \nActive cells: %i\nAll spatially tuned cells: %i (%.1f%%)\nPlace cells: %i (%.1f%%)\nLandmark cells: %i (%.1f%%)\nTrace cells: %i (%.1f%%)\n',sum(data.n_cells),sum(data.n_acs),sum([data.n_pcs,data.n_lcs,data.n_tcs]),(sum([data.n_pcs,data.n_lcs,data.n_tcs])/sum(data.n_acs))*100,sum(data.n_pcs),(sum(data.n_pcs)/sum(data.n_acs))*100,sum(data.n_lcs),(sum(data.n_lcs)/sum(data.n_acs))*100,sum(data.n_tcs),(sum(data.n_tcs)/sum(data.n_acs))*100);