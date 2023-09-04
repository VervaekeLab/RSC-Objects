function roiColoredByTag(roi_array, classes, varargin)
%% roiColoredByTag
% Plots the FOV of a session and color the ROIs within it according to the
% classes given as input.
%
% INPUT
%   -- Required
%   roi_array: A roi struct as provided from ROImanager containing information
%       about all rois.
%   classes: A int array like [0,0,1,2,0, ... ] where each number refers to
%       a class you want to color the rois. This must be of the same length
%       as the roi_array struct.
%
%   -- Optional
%   varargin: Alternating argument name and value, ex: ["keyword", value].
%       Arguments - See "Default parameters" section.
%
% Written by Andreas S Lande 2019

%% Default parameters
params.fov_image = []; % FOV image to be added in background.
params.display_method = "hardcoded"; % "plot" or "fill"

% Update parameters
params = sb.helper.updateParameterStruct(params,varargin);

%% Generate figure
figure(16);
clf;

n_classes = length(unique(classes));

% Find the number of classes
colormaps = cbrewer('qual','Set1',max(n_classes,3));

% If the fov image is not added, just use a mask of all rois
if isempty(params.fov_image)
    mask = zeros(size(roi_array(1).mask));
    imagesc(mask);
else
    imagesc(params.fov_image);
    colormap('gray');
end

hold on;

% For each ROI, draw the mask and colorcode based on class
for roi = 1:length(roi_array)
    
    boundary = bwboundaries(roi_array(roi).mask);
    boundary = boundary{1};
    
    if classes(roi) == 0
        roi_color = [0.5,0.5,0.5];
    else
        roi_color = colormaps(:,classes(roi));
    end
    
    switch params.display_method
        case "plot"
            plot(boundary(:,2),boundary(:,1),'Color',roi_color,'LineWidth',2);
        case "fill"
            if classes(roi) > 0
                f = fill(boundary(:,2),boundary(:,1),'k');
                f.FaceColor = roi_color;
                f.EdgeColor = roi_color;
                f.FaceAlpha = 1;
            end
            
        case "hardcoded"
            if classes(roi) > 0
                if classes(roi) == 1
                    roi_color = [239,35,60]./255;
                end
                
                if classes(roi) == 2
                    roi_color = [51,51,51]./255;
                end
                
                f = fill(boundary(:,2),boundary(:,1),'k');
                f.FaceColor = roi_color;
                f.EdgeColor = roi_color;
                f.FaceAlpha = 1;
                
            end
            
    end
    
end

xticks([]);
yticks([]);
axis equal





end
