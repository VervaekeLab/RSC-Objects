%% Load data
% Clean workspace
clear; 

% Close all open figures
close all;

% Load dataset
[mData,data] = sb.database.load.rsc_flyby();

% Parameters
signal_type = 'dff';

%% Analysis of only landmark cells
% Set variables
analysis = struct();
visualize = false;
counter = 1;

% Create analysis struct for each of the landmark cells
for s = 1:2:length(mData)-1
    
     % Find whisker side
    [~,w_ind] = max([length(mData(s+1).sData.landmarks.motor.right.on),length(mData(s+1).sData.landmarks.motor.left.on)]);
    if w_ind == 1
        side = 'right';
    else
        side = 'left';
    end
            
    % Get index for each roi that are landmark cells
    cell_inds = [mData(s).rmaps.lcs.ind,mData(s).rmaps.tcs.ind];
    
    for c = 1:length(cell_inds)
    
        % Get all rastermaps
        rmap_still = mData(s+1).rmaps.flyby.(signal_type).(side).still(:,:,cell_inds(c));
        rmap_moving = mData(s+1).rmaps.flyby.(signal_type).(side).moving(:,:,cell_inds(c));

        % Append index for session and roi
        analysis(counter).ind = cell_inds(c);
        analysis(counter).session = s;
        
        %Add rastermaps to struct
        analysis(counter).dff_motor_on = mData(s).rmaps.dff_motor_on(:,:,cell_inds(c));
        analysis(counter).dff_motor_off = mData(s).rmaps.dff_motor_off(:,:,cell_inds(c));
        analysis(counter).deconv_motor_on = mData(s).rmaps.deconv_motor_on(:,:,cell_inds(c));
        analysis(counter).deconv_motor_off = mData(s).rmaps.deconv_motor_off(:,:,cell_inds(c));
        
        % Run some analysis on the rastermaps for fly by
        bins_to_look_for_response = 135:155;
        response = sb.other.detectSignificantResponsesByShuffling(rmap_still(:,1:249),bins_to_look_for_response); % sample 31:185 is 2 seconds before event and 2 seconds after event 
        analysis(counter).still.response = response;
        analysis(counter).still.rmap = rmap_still;
        analysis(counter).still.rmap_deconv = mData(s+1).rmaps.flyby.deconv.(side).still(:,:,c);
        analysis(counter).still.rmap_dffZScored = mData(s+1).rmaps.flyby.dffZ.(side).still(:,:,c);
        
        response = sb.other.detectSignificantResponsesByShuffling(rmap_moving(:,1:249),bins_to_look_for_response); % sample 31:185 is 2 seconds before event and 2 seconds after event 
        analysis(counter).moving.response = response;
        analysis(counter).moving.rmap = rmap_moving;
        analysis(counter).moving.rmap_deconv = mData(s+1).rmaps.flyby.deconv.(side).moving(:,:,c);
        analysis(counter).moving.rmap_dffZScored = mData(s+1).rmaps.flyby.dffZ.(side).moving(:,:,c);

        analysis(counter).roiImage{1} = mData(s).sData.imdata.roiArray(cell_inds(c)).enhancedImage;
        analysis(counter).roiImage{2} = mData(s+1).sData.imdata.roiArray(cell_inds(c)).enhancedImage;
        
        % Plot the current roi
        if visualize
            
            % Find the max response across all rastermaps
            z_max = max([max(max(analysis(counter).dff_motor_on)),max(max(analysis(counter).dff_motor_off)),max(max(analysis(counter).still.rmap)),max(max(analysis(counter).moving.rmap))]);
    
            figure(1);clf;
            subplot(4,1,1);
            imagesc(analysis(counter).dff_motor_on,[0,z_max]);
            subplot(4,1,2);
            imagesc(analysis(counter).dff_motor_off,[0,z_max]);
            
            subplot(4,1,3);
            imagesc(analysis(counter).still.rmap,[0,z_max]);
            subplot(4,1,4);
            imagesc(analysis(counter).moving.rmap,[0,z_max]);
            colormap('jet');
            
        end
        
        % Increment counter
        counter = counter + 1;
    end
    
end

%% Analysis of all cells
% Create analysis_venn struct which contains data from all cells. This is
% used (among other things) to create the overall venn diagram.
counter = 1;
analysis_venn = struct();

for s = 1:2:length(mData)-1
    
    % Find whisker side
    [~,w_ind] = max([length(mData(s+1).sData.landmarks.motor.right.on),length(mData(s+1).sData.landmarks.motor.left.on)]);
    if w_ind == 1
        side = 'right';
    else
        side = 'left';
    end
    
    % Get index for each roi that are landmark cells
    cell_inds = 1:size(mData(s).rmaps.dff,3);
    
    
    for c = 1:length(cell_inds)
        
        % Get all rastermaps
        rmap_still = mData(s+1).rmaps.flyby.(signal_type).(side).still(:,:,c);
        rmap_moving = mData(s+1).rmaps.flyby.(signal_type).(side).moving(:,:,c);
        
        % Append index for session and roi
        analysis_venn(counter).ind = cell_inds(c);
        analysis_venn(counter).session = s;
        
        analysis_venn(counter).isLcs = ismember(c,unique([mData(s).rmaps.lcs.ind,mData(s).rmaps.tcs.ind]));
        analysis_venn(counter).isPcs = ismember(c,unique([mData(s).rmaps.pcs.ind]));
        
        
        %Add rastermaps to struct
        analysis_venn(counter).dff_motor_on = mData(s).rmaps.dff_motor_on(:,:,cell_inds(c));
        analysis_venn(counter).dff_motor_off = mData(s).rmaps.dff_motor_off(:,:,cell_inds(c));
        analysis_venn(counter).deconv_motor_on = mData(s).rmaps.deconv_motor_on(:,:,cell_inds(c));
        analysis_venn(counter).deconv_motor_off = mData(s).rmaps.deconv_motor_off(:,:,cell_inds(c));
        
        % Run some analysis on the rastermaps for fly by
        bins_to_look_for_response = 135:155; %110:125; % 31:185
        response = sb.other.detectSignificantResponsesByShuffling(rmap_still(:,1:249),bins_to_look_for_response); % sample 31:185 is 2 seconds before event and 2 seconds after event         analysis_venn(counter).still.response = response;
        analysis_venn(counter).still.response = response;
        analysis_venn(counter).still.rmap = rmap_still;
        analysis_venn(counter).still.rmap_deconv = mData(s+1).rmaps.flyby.deconv.(side).still(:,:,c);
        response = sb.other.detectSignificantResponsesByShuffling(rmap_moving(:,1:249),bins_to_look_for_response); % sample 31:185 is 2 seconds before event and 2 seconds after event
        
        analysis_venn(counter).moving.response = response;
        analysis_venn(counter).moving.rmap = rmap_moving;
        analysis_venn(counter).moving.rmap_deconv = mData(s+1).rmaps.flyby.deconv.(side).moving(:,:,c);
        
        analysis_venn(counter).roiImage{1} = mData(s).sData.imdata.roiArray(cell_inds(c)).enhancedImage;
        analysis_venn(counter).roiImage{2} = mData(s+1).sData.imdata.roiArray(cell_inds(c)).enhancedImage;
        
        % Increment counter
        counter = counter + 1;
    end
    
end

%% Summarize the number of responses
flyby_stats = struct();
flyby_stats.still = 0;
flyby_stats.moving = 0;
flyby_stats.both = 0;
flyby_stats.none = 0;

% For each cell, add a number to the correct category of response
for c = 1:length(analysis)
    
    if analysis(c).still.response
        if analysis(c).moving.response
            flyby_stats.both = flyby_stats.both + 1;
        else
            flyby_stats.still = flyby_stats.still + 1;
        end
        
    elseif analysis(c).moving.response
        flyby_stats.moving = flyby_stats.moving + 1;
    else
        flyby_stats.none = flyby_stats.none + 1;
    end
   
end

