function placeCells(sData,varargin)
%% placeCells
% Show all place cells found in the session.
%
% INPUT
%   -- Required
%       sData.
%
%   -- Optional
%   varargin: Alternating argument name and value, ex: ["keyword", value].
%       Arguments - See "Default parameters" section.
%
% Written by Andreas S Lande 2019


%% Default parameters
params = struct();

% Sort all rastermaps based on the motor manipulations
try
    params.manipulated_motor_position = sData.landmarks.motor.manipulated; % This can be overwritten as input
catch
    params.manipulated_motor_position = 0;
end

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


%% Generate data
pos = sData.behavior.wheelPosDsBinned;
lap = sData.behavior.wheelLapDsBinned;

if sum(sData.sessionInfo.sessionID(2:3) == '16') == 2
    dff = sData.imdata.roiSignals(2).dff;
else
    dff = sData.imdata.roiSignals(2).dffSubtractedNpil;
end

deconv = smoothdata(sData.imdata.roiSignals(2).deconv','gaussian',10)';

run = sData.behavior.runSpeedDs;

% Create dff rastermap
rmaps.dff = sb.generate.rasterMaps(pos,lap,dff);

% Create deconv rastermap
rmaps.deconv = sb.generate.rasterMaps(pos,lap,deconv,'smoothing',5);

% Create motor rastermap
right_motor = sData.daqdata.right_motor_speed(sData.daqdata.frameIndex);
rmaps.motor_right = sb.generate.rasterMaps(pos,lap,right_motor,'smoothing',0);

left_motor = sData.daqdata.left_motor_speed(sData.daqdata.frameIndex);
rmaps.motor_left = sb.generate.rasterMaps(pos,lap,left_motor,'smoothing',0);

% Create run speed map
rmaps.run = sb.generate.rasterMaps(pos,lap,run,'smoothing',0);

% If landmark is manipulated, sort rastermaps thereafter
if params.manipulated_motor_position
    rmaps.run = sb.sort.rasterMapByLandmarkOnset(rmaps.run,sData.landmarks.motor.on(:,params.manipulated_motor_position));
    rmaps.motor_right_sorted = sb.sort.rasterMapByLandmarkOnset(rmaps.motor_right,sData.landmarks.motor.on(:,params.manipulated_motor_position));
    rmaps.motor_left_sorted = sb.sort.rasterMapByLandmarkOnset(rmaps.motor_left,sData.landmarks.motor.on(:,params.manipulated_motor_position));
    rmaps.dff_sorted = sb.sort.rasterMapByLandmarkOnset(rmaps.dff,sData.landmarks.motor.on(:,params.manipulated_motor_position));
    rmaps.deconv_sorted = sb.sort.rasterMapByLandmarkOnset(rmaps.deconv,sData.landmarks.motor.on(:,params.manipulated_motor_position));
    
    % Find all place cells and sort their peaks
    pcs = sb.classify.placeTunedRoisSimple(rmaps.deconv_sorted(:,:,:));
    pcs = find(pcs);

    % Split rmaps into motor on vs motor off
    rmaps.dff_motor_on = rmaps.dff(find([sData.landmarks.motor.on(:,params.manipulated_motor_position) == 1]),:,:);
    rmaps.dff_motor_off = rmaps.dff(find([sData.landmarks.motor.on(:,params.manipulated_motor_position) == 0]),:,:);
    rmaps.deconv_motor_on = rmaps.deconv(find([sData.landmarks.motor.on(:,params.manipulated_motor_position) == 1]),:,:);
    rmaps.deconv_motor_off = rmaps.deconv(find([sData.landmarks.motor.on(:,params.manipulated_motor_position) == 0]),:,:);
  
    sb.plot.rasterMaps({rmaps.dff_motor_on(:,:,pcs),rmaps.dff_motor_off(:,:,pcs)},'add_lines_vertical',sData.landmarks.motor.pos,'smoothing',3);
    %sb.plot.rasterMaps({rmaps.deconv_motor_on(:,:,pcs),rmaps.deconv_motor_off(:,:,pcs)},'add_lines_vertical',sData.landmarks.motor.pos,'smoothing',3);
else
    
   
    pcs = sb.classify.placeTunedRoisSimple(rmaps.deconv(find([sData.landmarks.motor.on(:,1) == 1]),:,:));
    pcs = find(pcs);
    sb.plot.rasterMaps(rmaps.dff(:,:,pcs),'add_lines_vertical',[33,73]);
    
   
%     pcs = sb.classify.placeTunedRoisSimple(rmaps.deconv);
%     pcs = find(pcs);
%     sb.plot.rasterMaps(rmaps.dff(:,:,pcs));
%     
    
    
end




end