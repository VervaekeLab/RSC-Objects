%% This script creates an example figure of how the landmark manipulation
% PAPER 1 FIG 4
% is done and how it is sorted.

% Create a pseudo motor landmark map based on a session
s = 5;
landmark_width = 3;

landmarks = mData(5).sData.landmarks.motor.on(:,1);
landmarks = landmarks(find(~isnan(landmarks)));
sorted_landmarks = flipud(sort(landmarks));

% Set a value for each lap where a landmark is present
motormap = zeros(105,length(landmarks));
motormap(73:73+landmark_width-1,:) = ones(landmark_width,length(landmarks));
motormap_sorted = motormap;
motormap_sorted(33:33+landmark_width-1,1:sum(landmarks)) = ones(landmark_width,sum(landmarks));

for l = 1:length(landmarks)
   
    if landmarks(l) == 1
        motormap(33:33+landmark_width-1,l) = ones(1,landmark_width);
    end
    
end

% Plot
figure(1);
clf;

% Original
subplot(1,2,1);
imagesc(motormap');

cmap = colormap('gray');
colormap(flipud(cmap));
xticks([]);
yticks([]);


% Sorted
subplot(1,2,2);
imagesc(motormap_sorted');
xticks([]);
yticks([]);
hold on
line([1,105],[sum(landmarks),sum(landmarks)]+0.5,'Color','k');

set(gcf,'renderer', 'painters', 'Position', [2600 200 500 300])
