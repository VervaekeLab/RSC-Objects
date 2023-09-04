function sData = loadPupilData(sData)


sID = sData.sessionInfo.sessionID;

% Folder containing finished files
f_path = '/Users/annachristinagarvert/UIO Physiology Dropbox Dropbox/Lab Data/Andreas Lande/WHISKER CAM DATA/PupilVideo/PupilSignals';

% Load file
load(fullfile(f_path,[sID(1:17),'_pupilsignals.mat']));

seq_file_path = '/Users/annachristinagarvert/UIO Physiology Dropbox Dropbox/Lab Data/Andreas Lande/WHISKER CAM DATA/PupilVideo/';

% Load image stack to find first frame where the eye is bright
files_in_dir = dir(seq_file_path);
for f = 1:length(files_in_dir)
    if contains(files_in_dir(f).name,sID(1:17))
       seq_file = fullfile(seq_file_path,files_in_dir(f).name);
        break;
    end
end

IM = sb.database.load.getSelectedFrames(seq_file,1:1000);

first_frame = find(reshape(nanmean(IM(1,:,:),2),[1,1000])>10);
first_frame = first_frame(1);

% Include ondly the index where there is not nans in the radius signal,
% which is used as a proxy to know when there was actually imaging
include_bins = first_frame:size(Center,1);
pupil = struct(); 
pupil.x = Center(include_bins,1);
pupil.y = Center(include_bins,2);
pupil.height = Height(include_bins);
pupil.width = Width(include_bins);
pupil.orientation = Orientation(include_bins);
pupil.radius = Radius(include_bins);

sData.pupil = pupil;



end