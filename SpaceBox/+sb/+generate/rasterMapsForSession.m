function rmaps = rasterMapsForSession(sData, varargin)
% rasterMapsForSession
% Create a rastermap of all cells in sData
% INPUT
%   -- Required
%   sData
%
%   -- Optional
%   varargin: Alternating argument name and value, ex: ["keyword", value].
%       Arguments - See "Default parameters" section.
%
% OUTPUT
%   rmaps: Struct containing the rastermaps for dff, deconv, running speed,
%       all divided into motor on and motor off laps.


%% Default parameters
params = struct();
params.dff_signal = 'dffSubtractedNpil'; % Choose what type of dff signal to be used. Inside sData.imdata.roiSignals there is both dff and dffSubtractedNpil
params.run_threshold = 2; % bins where animal runs cm/s lower than this is excluded.
params.manipulated_motor = 1; % This means the first motor position is manipulated.
params.laps_on = [];
params.laps_off = [];

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

% Get signals
pos = sData.behavior.wheelPosDsBinned;
lap = sData.behavior.wheelLapDsBinned;
run = sData.behavior.runSpeedDs;
dff = sData.imdata.roiSignals(2).(params.dff_signal);
deconv = smoothdata(sData.imdata.roiSignals(2).deconv','gaussian',10)'; % Smooth the deconvolved signal for a 300 ms (-ish) window

% Create struct to store all rastermaps
rmaps = struct();

% Create a binary array of which sample is above or below a running speed threshold
run_thresholded = [run<params.run_threshold];
if sum(run_thresholded) > 10000; run_thresholded = zeros(size(run_thresholded));end

% Create rastermap for deconvolved signal
rmaps.deconv = sb.generate.rasterMaps(pos,lap,deconv,...
    'smoothing',5, ...
    'exclude_inds', run_thresholded);

% Create rastermap for dF/F signals
rmaps.dff = sb.generate.rasterMaps(pos,lap,dff, 'exclude_inds', run_thresholded);

% Create rastermap for run speed signal
rmaps.run = sb.generate.rasterMaps(pos,lap,run,'smoothing',0,'exclude_inds', run_thresholded);


% Create rastermap using the speed signal from the landmark motors
try
    right_motor = sData.daqdata.right_motor_speed(sData.daqdata.frameIndex);
    rmaps.motor_right = sb.generate.rasterMaps(pos,lap,right_motor,'smoothing',0);
    
    left_motor = sData.daqdata.left_motor_speed(sData.daqdata.frameIndex);
    rmaps.motor_left = sb.generate.rasterMaps(pos,lap,left_motor,'smoothing',0);
catch
    rmaps.motor_right = [];
    rmaps.motor_left = [];
end

% Divide rmaps into motor on and motor off laps if it has been manipulated
if params.manipulated_motor
    if isempty(params.laps_on)
        params.laps_on = find([sData.landmarks.motor.on(:,params.manipulated_motor) == 1]);
    end

    if isempty(params.laps_off)
        params.laps_off = find([sData.landmarks.motor.on(:,params.manipulated_motor) == 0]);
    end
else
    if isempty(params.laps_on)
        params.laps_on = 1:max(sData.behavior.wheelLapDs);
        params.laps_off = [];
    end
end

% Split the rastermaps into motor on and motor off.
rmaps.deconv_motor_on = rmaps.deconv(params.laps_on,:,:);
rmaps.deconv_motor_off = rmaps.deconv(params.laps_off,:,:);
rmaps.dff_motor_on = rmaps.dff(params.laps_on,:,:);
rmaps.dff_motor_off = rmaps.dff(params.laps_off,:,:);
rmaps.run_motor_on = rmaps.run(params.laps_on,:);
rmaps.run_motor_off = rmaps.run(params.laps_off,:);

% Add laps on and off to rmaps
rmaps.meta.laps_on = params.laps_on;
rmaps.meta.laps_off = params.laps_off;

end