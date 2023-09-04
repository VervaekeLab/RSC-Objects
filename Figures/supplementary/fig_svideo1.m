%% Clear
clear all;

%% Load images
whisker_img = stack2mat("/Volumes/Andreas1/Data/Paper_1_supp_video_raw/m1410_20191023_01_2PClockFrames/m1410_20191023_01_whiskerImages_02.tif"); % this loads image for frame 1001 until 2000.
whisker_img = cat(3,whisker_img,stack2mat("/Volumes/Andreas1/Data/Paper_1_supp_video_raw/m1410_20191023_01_2PClockFrames/m1410_20191023_01_whiskerImages_03.tif"));

calcium_img = stack2mat("/Volumes/Andreas1/Data/PROCESSED DATA/mouse1410/session-m1410-20191023-01_VR/calcium_images_aligned/calcium_images_m1410-20191023-01_VR_block001_plane001_ch2_part001.tif");
calcium_img = cat(3,calcium_img,stack2mat("/Volumes/Andreas1/Data/PROCESSED DATA/mouse1410/session-m1410-20191023-01_VR/calcium_images_aligned/calcium_images_m1410-20191023-01_VR_block001_plane001_ch2_part002.tif"));

% Smooth calcium images
ca_img_smoothed = zeros(size(calcium_img,1),size(calcium_img,2),size(calcium_img,3));
for x = 900:3100
    ca_img_smoothed(:,:,x) = nanmean(calcium_img(:,:,x-2:x+2),3);
end

calcium_img = ca_img_smoothed(:,:,1001:3000);

%% Load high speed whisker images
load("/Volumes/Andreas1/Data/Paper_1_supp_video_raw/supp_video_signals.mat");

whisker_img_high_fps = stack2mat("/Volumes/Andreas1/Data/Paper_1_supp_video_raw/supp_video_fast_frames.tif");

%% Combined the whisker images

skip = sum(signals.two_p_onset(1522:2200)==1);
whisker = cat(3,whisker_img(:,:,1:521),whisker_img_high_fps,whisker_img(:,:,521+skip+1:end));

% Frame rate signal
update_low_fps_data = signals.two_p_onset(1001:1000+size(whisker,3));

% 
% 
% whisker = whisker_img(:,:,1:500);
% c = 501;
% c_high = 1;
% for x = 1501:3001
%     
%     if signals.two_p_onset(x)
%         whisker = cat(3,whisker,whisker_img(:,:,c));
%         c = c+1;
%         
%     else
%        whisker = cat(3,whisker,whisker_img_high_fps(:,:,c_high));
%        c_high = c_high + 1;
%     end
%     
% end



%% Load session
data.sessionIDs = {'m1410_20191023_01'};
mData = sb.database.loadSessions(data.sessionIDs);

% Create rastermaps for each session
for s = 1:length(mData)
    mData(s).rmaps = sb.classify.cells(mData(s).sData,'manipulated_motor',1,'dff_signal','dffSubtractedNpil');
end

sData = mData(1).sData;
rmaps = mData(1).rmaps;
clearvars mData

%% Get data to plot
first_frame = 1000;
last_frame = first_frame + 31 * 60; % 31 fps for 60 seconds
pos = sData.behavior.wheelPosDsBinned(first_frame:end);
lap = sData.behavior.wheelLapDsBinned(first_frame:end);
lap = lap - lap(1) + 1;

% Rastermap
rmap = zScoreNormalize(rmaps.tcs.deconv_motor_on(:,:,6),'all');

% Selected cell
roi = 36;

% Create a signal to indicate the onset of landmark
landmark_sig = zeros(1,2000);

c = 1;
while c < 2000
    
    if ismember(pos(c),[32,33])
        landmark_sig(c) = 1;
        c = c+5;
    end
    
    if ismember(pos(c),[72,73])
        landmark_sig(c) = 1;
        c = c+5;
    end
    c = c+1;     
end

landmark_sig(landmark_sig<1) = -100;



%% Create a video of the calcium signal
figure(2); clf;
clearvars F

set(gcf,'renderer', 'painters', 'Position', [1800,1000,1030,1200]) %[400,1000,2400,610]

counter_low = 1;
counter_high = 1;

rmap = zeros(20,105);
frames = first_frame:last_frame;
    
for low_fps_onset = update_low_fps_data
    
    figure(2);
%      col = 2;
%      row = 5;
%     for xx = 1:col*row
%        subplot(row,col,xx);
%        imagesc(ones(10,10));
%        text(5,5,sprintf('%i',xx),'FontSize',14);
%     end
    
    
    % Clear figure
    if low_fps_onset
        clf;
    end
    
    
    % Video of whiskers
    subplot(row,col,[1:4])
    imshow(squeeze(imresize(whisker(:,:,counter_high),472/400))); % offset by 10 to match the first frame at 1011. 
    colormap('gray');
    hold on;
    
    if ismember(counter_low,[521:521+681])
        title('Slow motion playback (300 fps)','FontSize',18,'Color','r');
    else
        title('Real-time playback (31 fps)','FontSize',18,'Color','k');
    end

    counter_high = counter_high+1;
    
    % Set color
    ax = subplot(row,col,[6,8]);
    colormap(ax,[1,1,1;0,0.1,1])
    
    if low_fps_onset
        
        current_low_fps_signal_frame = 1000+counter_low;

        % Two photon FOV
        subplot(row,col,[5,7]);

        % Create a FOV of the imaging where the roi is superimposed as a small image
        small_roi_fov = calcium_img(343:382,350:389,counter_low);
        zoom_factor = 4;
        small_roi_fov = imresize(small_roi_fov,zoom_factor);

        fov_to_plot = calcium_img(:,:,counter_low);

        offset = 10;
        fov_to_plot(offset:offset+size(small_roi_fov,1)-1,offset:offset+size(small_roi_fov,2)-1) = small_roi_fov;

        % Add white border
        fov_to_plot(offset:offset+size(small_roi_fov,1)-1,[offset-1:offset]) = 255;
        fov_to_plot(offset:offset+size(small_roi_fov,1)-1,[size(small_roi_fov,1)+offset-1:size(small_roi_fov,1)+offset]) = 255;
        fov_to_plot([offset-1:offset],offset:offset+size(small_roi_fov,1)-1) = 255;
        fov_to_plot([size(small_roi_fov,1)+offset-1:size(small_roi_fov,1)+offset],offset:offset+size(small_roi_fov,1)-1) = 255;

        imshow(fov_to_plot,[0,180]);
        hold on; 
        viscircles([23*zoom_factor,23*zoom_factor],6*zoom_factor,'EnhanceVisibility',false,'Color','w','LineWidth',1.4);
        viscircles(sData.imdata.roiArray(36).center,6,'EnhanceVisibility',false,'Color','w','LineWidth',1.4);
        
        % dF/F signal
        subplot(row,col,[9,10])
        plot(sData.imdata.roiSignals(2).dffSubtractedNpil(roi,first_frame:current_low_fps_signal_frame),'Color','k'); 
        hold on;
        deconv_signal = sData.imdata.roiSignals(2).deconv(roi,first_frame:current_low_fps_signal_frame);
        deconv_signal(deconv_signal>0) = 0.5;
        deconv_signal(deconv_signal<0.4) = nan;
        plot(deconv_signal,'Color','b','Marker','.','MarkerSize',18); 
        line([counter_low,counter_low],[-0.2,1.5],'Color',[0.6,0.6,0.6]);
        plot(landmark_sig,'Color','r');
        hold off;
        ylim([-0.2,1])
        xlim([1,last_frame-first_frame])
        xticks([1,31*20,31*40,31*60]);
        xticklabels({0,20,40,60})
        yticks([]);
        xlabel('Time (seconds)');
        set(gca,'FontSize',14)

        % Raster map of the cells response
        ax = subplot(row,col,[6,8]);

        if rmap(lap(counter_low),pos(counter_low)) == 0
            if deconv_signal(counter_low)>0
                rmap(lap(counter_low),pos(counter_low)) = 1;
            end
        end

        imagesc(rmap,[0,1]);
        colormap(ax,[1,1,1;0,0.1,1])
        xlim([0,105.5])
        hold on;
        line([33,33],[0,51],'LineWidth',2,'Color','r')
        line([73,73],[0,51],'LineWidth',2,'Color','r')
        line([93,93],[0,51],'LineWidth',2,'Color','k')
        plot(pos(counter_low),lap(counter_low),'LineStyle','none','Marker','.','MarkerSize',20,'Color','k')
        ylim([0.5,18.5])
        xticks([0,105])
        xticklabels({0,157})
        ylabel('Lap #');
        xlabel('Position (cm)');
        %title(sprintf('CL: %i, CH: %i',counter_low,counter_high));
        
        % Create average of rmap
        avg_sig = normalize(smoothdata(nanmean(rmap),'gaussian',5),'range',[0,3]);

        plot(18.5-avg_sig,'Color','k','LineWidth',2);

        set(gca,'FontSize',14)
        counter_low = counter_low+1;
    end
    
    % Capture frame
    F(counter_high) = getframe(gcf);
    drawnow
    

end

F = F(2:end);

% create the video writer with 1 fps
writerObj = VideoWriter('signalVideo.avi');
writerObj.FrameRate = 31;

% open the video writer
open(writerObj);
% write the frames to the video
for current_low_fps_signal_frame=1:length(F)
    % convert the image to a frame
    frame = F(current_low_fps_signal_frame) ;
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);


















