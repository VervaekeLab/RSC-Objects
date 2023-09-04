

% these are hard coded, the landmarks ID in lcsA, lcsB can be variable due
% to the randomness the shuffling causes,these cells were chosen from
% cells that were identified landmark cells, find the classification here:
% /Users/annachristinagarvert/UIO Physiology Dropbox Dropbox/Lab Data/Ann Christin Garvert/Andreas/MichelesData/CueTuned_shuffle

% sanpaper 
% session 3 cellIdx: 20   --> 'm5115-20201202-01' :122
% session 9 cellIdx: 3, 8 --> 'm5108-20201120-01' : 66, 242

% mixed in sandpaper
% session 3 cellIdx: 4 --> 'm5115-20201202-01' :13

% felt pad
% session 2: cellIdx = 2,5,12,16; --> 'm5115-20201203-01': 6, 11, 29, 45
% session 4: cellIdx: 11,14 ,15   --> 'm5114-20201203-01': 51, 14,  15
% mixed in feltpads = session 3, celldIdx : 20  --> 'm5115-20201202-01' :84

%% Plot example place cell, landmark cell and trace cell
% All from session 1.

% Landmark cells
% Plot two cells that are felt pad cells
plotCells(cat(3,mData(4).rmaps.deconv(1:60,:,[25,51]),mData(2).rmaps.deconv(1:60,:,6),...
     mData(5).rmaps.deconv(1:60,:,79)))

plotCells(cat(3,mData(9).rmaps.deconv(1:60,:,[66,127]),...
    mData(3).rmaps.deconv(1:60,:,[32,122])))
% mData(3).rmaps.deconv(1:60,:,[16,3]),mData(6).rmaps.deconv(1:60,:,[53]) mData(9).rmaps.deconv(1:60,:,[313])

plotCells(cat(3, mData(3).rmaps.deconv(1:60,:,[3 97 103 191])))

function plotCells(rmaps)

    % Parameters
    z_score_level = 4;

    figure;
    clf;
    landmark_cells_color_1 = [101,191,164]/255;
    landmark_cells_color_2 = [140,160,196]/255;

    for r = 1:size(rmaps,3)
        subplot(size(rmaps,3),1,r)
        
        rmap_to_plot = normalize(rmaps(10:60,:,r),2);
        %         rmap_to_plot = zScoreNormalize(rmaps(10:40,:,r),'all');

        imagesc(rmap_to_plot,[0,z_score_level]);
        colormap(flipud(colormap('gray')));

        % Draw landmark lines
        hold on
        
        cue_1 = [42 47 122 128]./1.5;
        cue_2 = [63 68 100 106]./1.5;

        for i = 1:length(cue_1)
                line([cue_1(i),cue_1(i)],[0,100],'LineWidth',0.5,'Color',landmark_cells_color_1, 'LineStyle','-','LineWidth',1);
                line([cue_2(i),cue_2(i)],[0,100],'LineWidth',0.5,'Color',landmark_cells_color_2, 'LineStyle','-','LineWidth',1);
        end
        %line([93,93],[0,100],'LineWidth',1,'Color',[0.7,0.7,0.7]);
        
        yticks([1,50])
        if r == size(rmaps,3)
            xticks([1,105]);
            xticklabels({0,157})
        else
            xticks([]);
        end

        set(gca,'FontSize',22)
    end

    set(gcf,'Renderer', 'painters', 'Position', [1200 200 350 900])
end