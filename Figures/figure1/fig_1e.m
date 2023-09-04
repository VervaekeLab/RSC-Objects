%% Plot example place cell, landmark cell and trace cell
% All from session 1.
session = 1;

% Landmark cells
% Plot two cells that are trace landmark cells and two that are regular
% landmark cells
plotCells(cat(3,mData(session).rmaps.tcs.deconv_motor_on(:,:,[13,17]),mData(session).rmaps.lcs.deconv_motor_on(:,:,[8,16])))

% Place cell
i = [49,9,57,39]; % index of 4 place cells
plotCells(mData(session).rmaps.pcs.deconv_motor_on(:,:,i));

function plotCells(rmaps)

    % Parameters
    z_score_level = 4;

    figure;
    clf;
    landmark_cells_color = [239,35,60]/255;
    

    for r = 1:size(rmaps,3)
        subplot(size(rmaps,3),1,r)
        
        rmap_to_plot = normalize(rmaps(10:40,:,r),2);
        %         rmap_to_plot = zScoreNormalize(rmaps(10:40,:,r),'all');

        imagesc(rmap_to_plot,[0,z_score_level]);
        colormap(flipud(colormap('gray')));

        % Draw landmark lines
        hold on
        line([32,32],[0,100],'LineWidth',0.5,'Color',landmark_cells_color, 'LineStyle','-','LineWidth',1);
        line([72,72],[0,100],'LineWidth',0.5,'Color',landmark_cells_color, 'LineStyle','-','LineWidth',1);
        %line([93,93],[0,100],'LineWidth',1,'Color',[0.7,0.7,0.7]);
        
        yticks([1,30])
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