function flybyResponses(sData,varargin)
%% flybyResponses
%
% INPUT
%   -- Required
%      sData
%
%   -- Optional
%   varargin: Alternating argument name and value, ex: ["keyword", value].
%       Arguments - See "Default parameters" section.
%
% OUTPUT
%
%
% Written by Andreas S Lande 2021


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

% Update parameters
params = sb.helper.updateParameterStruct(params,varargin);



% Get motor onsets
run = sData.behavior.runSpeedDs;
left_motor_sig = sData.daqdata.left_motor_speed(sData.daqdata.frameIndex);
right_motor_sig = sData.daqdata.right_motor_speed(sData.daqdata.frameIndex);

[right_onsets_still,right_onsets_moving] = sb.other.findLandmarkFlybyOnsets(run,right_motor_sig);
[left_onsets_still,left_onsets_moving] = sb.other.findLandmarkFlybyOnsets(run,left_motor_sig);


% Create rastermap of responses for each cell in a time window around the
% touch presentation

dff = sData.imdata.roiSignals(2).dffSubtractedNpil;
window = [-(31*4):(31*5)];

% Create rastermap for right motor
rmaps_right = zeros(length(right_onsets_still),length(window),size(dff,1));

for c = 1:size(dff,1)
    for trial = 1:length(right_onsets_still)
        rmaps_right(trial,:,c) = dff(c,right_onsets_still(trial) + window);
    end
end

% Create rastermap for left motor
rmaps_left = zeros(length(left_onsets_still),length(window),size(dff,1));

for c = 1:size(dff,1)
    for trial = 1:length(left_onsets_still)
        rmaps_left(trial,:,c) = dff(c,left_onsets_still(trial) + window);
    end
end

close all;
sb.plot.rasterMaps({rmaps_right,rmaps_left},'smoothing',3,'add_lines_vertical',124);

