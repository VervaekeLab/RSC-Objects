function [landmarkOnsets,motorSignal] = findMotorOnOffLapsFromWhiskerImaging(sData,varargin)

% Variables
pos = sData.behavior.wheelPosDsBinned;
lap = sData.behavior.wheelLapDsBinned;

if length(sData.whiskers.rm_frames)>length(sData.whiskers.lm_frames)
    % The rm / lm notation is reversed, so rm is left and lm is right...
    motor = sData.whiskers.left_motor_sig;
else
    motor = sData.whiskers.right_motor_sig;
end

motor = normalize(-motor,'range',[0,100]);
motor(motor<40) = 0;
motor(motor>90) = 0;
motor = smoothdata(motor,'movmean',20);

[~,LEDOnsets] = findpeaks(sData.whiskers.led_sig,'MinPeakDistance',4,'MinPeakHeight',200,'MinPeakWidth',1);

motorSignal = zeros(1,length(pos));
motorSignal(1:length(LEDOnsets)) = motor(LEDOnsets);
fprintf('\nFound %i samples where LED is on. Total imaged frames is %i\n', length(LEDOnsets),length(pos));

motorRasterMap = sb.generate.rasterMaps(pos,lap,motorSignal);
%sb.plot.rasterMaps(motorRasterMap);

first_landmark = nanmean(motorRasterMap(:,20:40)');
first_landmark(first_landmark>0) = 1;

second_landmark = nanmean(motorRasterMap(:,60:80)');
second_landmark(second_landmark>0) = 1;

landmarkOnsets = [];
landmarkOnsets(:,1) = first_landmark;
landmarkOnsets(:,2) = second_landmark;

end


