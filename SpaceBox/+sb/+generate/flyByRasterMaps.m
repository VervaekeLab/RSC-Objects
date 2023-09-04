function rmaps = flyByRasterMaps(sData,varargin)
%% flyByRasterMaps
% Creates rastermaps for a flyby experiment
% INPUT
%   -- Required
%   sData:
%
%   -- Optional
%   varargin: Alternating argument name and value, ex: ["keyword", value].
%       Arguments - See "Default parameters" section.
%
% OUTPUT
%   


%% Default parameters
params = struct();
params.signal_type = 'dff';

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

% Find flyby onsets
[right_onsets_still,right_onsets_moving] = sb.other.findLandmarkFlybyOnsets(run,right_motor_sig);
[left_onsets_still,left_onsets_moving] = sb.other.findLandmarkFlybyOnsets(run,left_motor_sig);

% Create rastermap of responses for each cell in a time window around the
% touch presentation
switch params.signal_type
    case 'dff'
        dff = sData.imdata.roiSignals(2).dffSubtractedNpil;
    case 'deconv'
        dff = sData.imdata.roiSignals(2).deconv;
        
    case 'dffZScored'
       dff = sData.imdata.roiSignals(2).dffZScored;
end

window = [-(31*4):(31*5)];

% ---- Create rastermaps for STILL laps
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

rmaps.left.still = rmaps_left;
rmaps.right.still = rmaps_right;

try
% ---- Create rastermaps for MOVING laps
% Create rastermap for right motor
rmaps_right = zeros(length(right_onsets_moving),length(window),size(dff,1));

for c = 1:size(dff,1)
    for trial = 1:length(right_onsets_moving)
        rmaps_right(trial,:,c) = dff(c,right_onsets_moving(trial) + window);
    end
end

% Create rastermap for left motor
rmaps_left = zeros(length(left_onsets_moving),length(window),size(dff,1));

for c = 1:size(dff,1)
    for trial = 1:length(left_onsets_moving)
        rmaps_left(trial,:,c) = dff(c,left_onsets_moving(trial) + window);
    end
end

rmaps.left.moving = rmaps_left;
rmaps.right.moving = rmaps_right;

catch

    rmaps.left.moving = [];
    rmaps.right.moving = [];

end


end

