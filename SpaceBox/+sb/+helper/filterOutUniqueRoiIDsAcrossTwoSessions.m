function [sData1, sData2] = filterOutUniqueRoiIDsAcrossTwoSessions(sData1,sData2,varargin)
%% filterOutUniqueRoiIDsAcrossTwoSessions
% This functions returns the sData structs given as input but after
% removing all rois that don't have a unique ID present in both sessions.
%
% INPUT
%   -- Required
%   sData1: sData from a session
%   sData2: sData from another session to compare with sData1.
%
%   -- Optional
%   varargin: Alternating argument name and value, ex: ["keyword", value].
%       Arguments - See "Default parameters" section.
%
% OUTPUT
%   sData1: sData from a session, where only rois also present in sData2 is
%       included.
%   sData2: sData from another session, where only rois also present in 
%       sData1 is included.


%% Default parameters
params = struct();

% Print all variables if "help" is only input 
try
    if inp == "help"
        global verbose; verbose = true;
        sb.helper.displayParameters(params)
        verbose = false;
        return;
    end
end


% Sort so that sData1 is the sData with most rois? 
% pass

% Update parameters
params = sb.helper.updateParameterStruct(params,varargin);

uIDs1 = {};
for r = 1:length(sData1.imdata.roiArray)
   uIDs1{end+1} = sData1.imdata.roiArray(r).uid; 
end

uIDs2 = {};
for r = 1:length(sData2.imdata.roiArray)
   uIDs2{end+1} = sData2.imdata.roiArray(r).uid; 
end

rois_to_include_1 = ismember(uIDs1,uIDs2);
rois_to_include_2 = ismember(uIDs2,uIDs1);

% Update roiSignals
all_fields = fields(sData1.imdata.roiSignals(2));
for f = 1:length(all_fields)
    if ~strcmp(all_fields{f},'ch')
    	sData1.imdata.roiSignals(2).(all_fields{f}) = sData1.imdata.roiSignals(2).(all_fields{f})(rois_to_include_1,:);
    end
end

% Update roiSignals
all_fields = fields(sData2.imdata.roiSignals(2));
for f = 1:length(all_fields)
    if ~strcmp(all_fields{f},'ch')
    	sData2.imdata.roiSignals(2).(all_fields{f}) = sData2.imdata.roiSignals(2).(all_fields{f})(rois_to_include_2,:);
    end
end


% Update roiArray
sData1.imdata.roiArray = sData1.imdata.roiArray(rois_to_include_1);
sData2.imdata.roiArray = sData2.imdata.roiArray(rois_to_include_2);

% Verify that it works
agreed = zeros(1,length(sData1.imdata.roiArray));

for r = 1:length(sData1.imdata.roiArray)
   
    uid = sData1.imdata.roiArray(r).uid;
    
    % CHeck that Uid is in sData2
    
    for r2 = 1:length(sData2.imdata.roiArray)
       uid2 = sData2.imdata.roiArray(r2).uid;
       
       if strcmp(uid,uid2)
          agreed(r) = 1; 
          break;
       end
        
    end
    
end

if length(sData1.imdata.roiArray) == sum(agreed)
    disp('ROI arrays across two sessions have been filtered');
else
   error('ROI arrays do not match') 
end








