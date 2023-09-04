


clear; close all;

%Requires dataset rsc_manipulation
[mData,data] = sb.database.load.rsc_manipulation;
    
firstSession = [1,3:6,8,10,13,14,18];
secondSession= [1,3,4,8,9,12,19,25,29,31];
fourthSession= [1,2,5,11,12,14,16,18,21,23];
plotCells(mData(1).rmaps.lcs.deconv_motor_on(:,:,firstSession));
plotCells(mData(2).rmaps.lcs.deconv_motor_on(:,:,secondSession));
plotCells(mData(4).rmaps.lcs.deconv_motor_on(:,:,fourthSession));



function plotCells(rmaps)

    % Parameters
    z_score_level = 4;

    figure;
    clf;
    landmark_cells_color = [239,35,60]/255;
    

    for r = 1:size(rmaps,3)
        subplot(10,1,r)
        
        rmap_to_plot = normalize(rmaps(10:60,:,r),2);
        %         rmap_to_plot = zScoreNormalize(rmaps(10:40,:,r),'all');

        imagesc(rmap_to_plot,[0,z_score_level]);
        colormap(flipud(colormap('gray')));

        % Draw landmark lines
        hold on
        line([32,32],[0,100],'LineWidth',0.5,'Color',landmark_cells_color, 'LineStyle','-','LineWidth',1);
        line([72,72],[0,100],'LineWidth',0.5,'Color',landmark_cells_color, 'LineStyle','-','LineWidth',1);
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