function sData = convertWhiskerImagingToLandmarkOnsetVector(sData)

params.visualize = false;

pos = sData.behavior.wheelPosDsBinned;
lap = sData.behavior.wheelLapDsBinned;

if length(sData.whiskers.rm_frames) > length(sData.whiskers.lm_frames)
    whisker_im_frames = [sData.whiskers.rm_frames(:).frame_num_imaging];
else
    whisker_im_frames = [sData.whiskers.lm_frames(:).frame_num_imaging];
end
whisker_im_frames = whisker_im_frames(1:4:end);

whisker_signal = zeros(1,length(lap));
whisker_signal(whisker_im_frames) = 1;

whiskmap = sb.generate.rasterMaps(pos,lap,whisker_signal);

first_motor = nanmean(whiskmap(:,30:45)');
    first_motor(first_motor>0) = 1;
    
if params.visualize
    sb.plot.rasterMaps(whiskmap);
    figure;
    plot(first_motor)
end

right_motor = [first_motor;ones(1,length(first_motor))];
sData.landmarks.motor.on = right_motor';


end