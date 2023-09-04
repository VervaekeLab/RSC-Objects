function allCells(sData,varargin)
%% allCells
% Plot all cells "as is" for a session
%
% INPUT
%   -- Required
%       sData
%
%   -- Optional
%   varargin: Alternating argument name and value, ex: ["keyword", value].
%       Arguments - See "Default parameters" section.
%
% OUTPUT
%
%

%% Default parameters
params = struct();
params.signal = "dffSubtractedNpil";
params.sort_by_motor_presence = 0;

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

%% Create rastermaps

pos = sData.behavior.wheelPosDsBinned;
lap = sData.behavior.wheelLapDsBinned;
run = sData.behavior.runSpeedDs;

signal = sData.imdata.roiSignals(2).(params.signal);

rmaps = sb.generate.rasterMaps(pos,lap,signal,'exclude_inds',run<2);

if params.sort_by_motor_presence
    
    
    motor_on_laps = find(sData.landmarks.motor.on(:,params.sort_by_motor_presence) == 1);
    motor_off_laps = find(sData.landmarks.motor.on(:,params.sort_by_motor_presence) == 0);
    
    rmaps_motor_on = rmaps(motor_on_laps,:,:);
    rmaps_motor_off = rmaps(motor_off_laps,:,:);
    
    % Plot
    sb.plot.rasterMaps({rmaps_motor_off,rmaps_motor_on},'add_lines_vertical',[31,71,91],'smoothing',1);
else
    % Plot
    sb.plot.rasterMaps(rmaps,'add_lines_vertical',[31,71,91]); 
end










end