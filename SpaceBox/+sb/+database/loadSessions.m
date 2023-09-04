function mData = loadSessions(sessionIDs)
%% loadSessions
%
% INPUT
%   -- Required
%       sessionIDs: Cell array where each element is the sessionID of a
%           session to be loaded. 
%           Ex: {'m1410_20191023_01','m1424_20200311_01'}.
%
% OUTPUT
%   mData: A structure containing the sData for each recording.


%% Load all sessions
mData = struct();

% Add .mat to end of all sessionIDs
for f = 1:length(sessionIDs)
   sessionIDs{f} = [sessionIDs{f},'.mat'];
    
end

% Load each sData
for f = 1:length(sessionIDs)
    load(fullfile(getPathToDir('experimentdata'),sessionIDs{f}));
    mData(f).sData = sData;
end