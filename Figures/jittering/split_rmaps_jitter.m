function [rmaps,motors] = split_rmaps_jitter(sData, varargin)
% rasterMapsForSession
% Create a rastermap of all cells in sData
% INPUT
%  -- Required
%  sData
%
%  -- Optional
%  varargin: Alternating argument name and value, ex: [“keyword”, value].
%    Arguments - See “Default parameters” section.
%
% OUTPUT
%  rmaps: Struct containing the rastermaps for dff, deconv, running speed,
%    all divided into motor on and motor off laps.

%% Default parameters
params = struct();
params.run_threshold = 2; % bins where animal runs cm/s lower than this is excluded.
params.manipulated_motor = 1; % This means the first motor position is manipulated.

motors = struct();
motors.laps_on = [];
motors.laps_off = [];
motors.ind.laps_off = [];

motors.laps_firstpair = [];
motors.ind.firstpair = [];

motors.laps_secondpair = [];
motors.ind.secondpair = [];

motors.laps_thirdpair = [];
motors.ind.thirdpair = [];

motors.laps_fourthpair = [];
motors.ind.fourthpair = [];

% Print all variables if “help” is only input
% try
%   if inp == “help”
%     global verbose; verbose = true;
%     sb.helper.displayParameters(params)
%     verbose = false;
%     return;
%   end
% end

% Update parameters
params = sb.helper.updateParameterStruct(params,varargin);

% Get signals
pos = sData.behavior.wheelPosDsBinned;
lap = sData.behavior.wheelLapDsBinned;
run = sData.behavior.runSpeedDs;
dff = sData.imdata.roiSignals(2).dff; 
deconv = smoothdata( sData.imdata.roiSignals(2).deconv, 'gaussian', 10);

% Create struct to store all rastermaps
rmaps = struct();

% Create a binary array of which sample is above or below a running speed threshold
run_thresholded = [run<params.run_threshold];

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

%Set bin for each pair
bin_firstpair_firstloc = 22;
bin_firstpair_secondloc = 55;
bin_secondpair_firstloc = 32;
bin_secondpair_secondloc = 65;
bin_thirdpair_firstloc = 41;
bin_thirdpair_secondloc = 75;
bin_fourthpair_firstloc = 52;
bin_fourthpair_secondloc = 85;

%Remove NaN
rmaps.motor_right(isnan(rmaps.motor_right)) = 0;

%Divide rmaps into laps based on motor location (for jittering)
if params.manipulated_motor
    for lap = 1:size(rmaps.motor_right,1)
      if rmaps.motor_right(lap,bin_firstpair_firstloc) && rmaps.motor_right(lap, bin_firstpair_secondloc) > 0
        if isempty(motors.laps_firstpair)
          motors.laps_firstpair = rmaps.motor_right(lap,:); 
        else
          motors.laps_firstpair = cat(1, motors.laps_firstpair, rmaps.motor_right(lap,:));
        end
        motors.ind.firstpair = cat(1, motors.ind.firstpair, lap);
      elseif rmaps.motor_right(lap,bin_secondpair_firstloc) && rmaps.motor_right(lap, bin_secondpair_secondloc) > 0
        if isempty(motors.laps_secondpair)
          motors.laps_secondpair = rmaps.motor_right(lap,:);
        else
          motors.laps_secondpair = cat(1,motors.laps_secondpair, rmaps.motor_right(lap,:));
        end
        motors.ind.secondpair = cat(1, motors.ind.secondpair, lap);
      elseif rmaps.motor_right(lap,bin_thirdpair_firstloc) && rmaps.motor_right(lap, bin_thirdpair_secondloc) > 0
        if isempty(motors.laps_thirdpair)
          motors.laps_thirdpair = rmaps.motor_right(lap,:);
        else
          motors.laps_thirdpair = cat(1, motors.laps_thirdpair, rmaps.motor_right(lap,:));
        end
        motors.ind.thirdpair = cat(1, motors.ind.thirdpair, lap);
      elseif rmaps.motor_right(lap,bin_fourthpair_firstloc) && rmaps.motor_right(lap, bin_fourthpair_secondloc) > 0
        if isempty(motors.laps_fourthpair)
          motors.laps_fourthpair = rmaps.motor_right(lap,:);
        else
          motors.laps_fourthpair = cat(1, motors.laps_fourthpair, rmaps.motor_right(lap,:));
        end
        motors.ind.fourthpair = cat(1, motors.ind.fourthpair,lap);
      else
        if isempty(motors.laps_off)
          motors.laps_off = rmaps.motor_right(lap,:);
        else
          motors.laps_off = cat(1, motors.laps_off, rmaps.motor_right(lap, :));
        end
        motors.ind.laps_off = cat(1, motors.ind.laps_off, lap);
      end
    end
end

% Divide rmaps into motor on and motor off laps if it has been manipulated
%if params.manipulated_motor
%  if isempty(params.laps_on)
%    params.laps_on = find([sData.landmarks.motor.on(:,params.manipulated_motor) == 1]);
%  end
%
%  if isempty(params.laps_off)
%    params.laps_off = find([sData.landmarks.motor.on(:,params.manipulated_motor) == 0]);
%  end
%else
%  if isempty(params.laps_on)
%    params.laps_on = 1:max(sData.behavior.wheelLapDs);
%    params.laps_off = [];
%  end
%end

% Split the rastermaps into motor on and motor off.
%rmaps.deconv_motor_on = rmaps.deconv(params.laps_on,:,:);
%rmaps.deconv_motor_off = rmaps.deconv(params.laps_off,:,:);
%rmaps.dff_motor_on = rmaps.dff(params.laps_on,:,:);
%rmaps.dff_motor_off = rmaps.dff(params.laps_off,:,:);
%rmaps.run_motor_on = rmaps.run(params.laps_on,:);
%rmaps.run_motor_off = rmaps.run(params.laps_off,:);

% Split the rastermaps based on location of motors
rmaps.deconv_firstpair = rmaps.deconv(motors.ind.firstpair,:,:);
rmaps.deconv_secondpair = rmaps.deconv(motors.ind.secondpair, :, :);
rmaps.deconv_thirdpair = rmaps.deconv(motors.ind.thirdpair, :, :);
rmaps.deconv_fourthpair = rmaps.deconv(motors.ind.fourthpair, :, :);
rmaps.deconv_motor_off = rmaps.deconv(motors.ind.laps_off,:,:);

rmaps.dff_firstpair = rmaps.dff(motors.ind.firstpair,:,:);
rmaps.dff_secondpair = rmaps.dff(motors.ind.secondpair, :, :);
rmaps.dff_thirdpair = rmaps.dff(motors.ind.thirdpair, :, :);
rmaps.dff_fourthpair = rmaps.dff(motors.ind.fourthpair, :, :);
rmaps.dff_motor_off = rmaps.dff(motors.ind.laps_off,:,:);

%Split rastermaps into motor on and motor off 
%rmaps.dff_motor_on = rmaps.dff(31:size(rmaps.motor_right,1), :,:);
%rmaps.deconv_motor_on = rmaps.deconv(31:size(rmaps.motor_right,1), :,:);

rmaps.dff_motor_off = rmaps.dff(motors.ind.laps_off,:,:);
rmaps.deconv_motor_off = rmaps.deconv(motors.ind.laps_off,:,:);

%!!NB if motors were introduced from first lap, remember to use these lines instead
rmaps.dff_motor_on = rmaps.dff(:, :, :);
rmaps.deconv_motor_on = rmaps.deconv(:, :, :);

% Add laps on and off to rmaps
%rmaps.meta.laps_on = params.laps_on;
rmaps.meta.laps_off = motors.laps_off;

end