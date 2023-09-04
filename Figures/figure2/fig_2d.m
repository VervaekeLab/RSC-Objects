%% Find all the cells that are within each group
only_engaged = [];
only_still = []; % this is renamed to be passive stimulus in the paper
only_moving = []; % This is renamed to "unengaged" in the paper
engaged_rmaps = [];
still_rmaps = [];
moving_rmaps = [];
signal_type = 'dff_motor_on';

% Only engaged
for i = 1:length(analysis_venn)
    
    % Engaged cells
    if analysis_venn(i).isLcs
        
        if sum([analysis_venn(i).still.response, analysis_venn(i).moving.response]) == 0

            if isempty(engaged_rmaps)
                engaged_rmaps = sb.generate.averagedTuningFromRasterMaps(analysis_venn(i).(signal_type));
            else
                engaged_rmaps = cat(1, engaged_rmaps, sb.generate.averagedTuningFromRasterMaps(analysis_venn(i).(signal_type)));
            end
        end
    else
        
        % Only passive (immobile)
        if [analysis_venn(i).still.response == 1] && [analysis_venn(i).moving.response == 0]
 
            if isempty(still_rmaps)
                still_rmaps = sb.generate.averagedTuningFromRasterMaps(analysis_venn(i).(signal_type));
            else
                still_rmaps = cat(1, still_rmaps, sb.generate.averagedTuningFromRasterMaps(analysis_venn(i).(signal_type)));
            end

        end
        
        % Only passive (running)
        if [analysis_venn(i).still.response == 0] && [analysis_venn(i).moving.response == 1]
            
            if isempty(moving_rmaps)
                moving_rmaps = sb.generate.averagedTuningFromRasterMaps(analysis_venn(i).(signal_type));
            else
                moving_rmaps = cat(1, moving_rmaps, sb.generate.averagedTuningFromRasterMaps(analysis_venn(i).(signal_type)));
            end
            
        end
        
    end
 
end



%% Create bar plot
figure(56); clf;

lm1 = 32:43;
lm2 = 72:83;

engaged_max_dff_at_landmark = nanmean(engaged_rmaps(:,[lm1,lm2])');
still_max_dff_at_landmark = nanmean(still_rmaps(:,[lm1,lm2])');
moving_max_dff_at_landmark = nanmean(moving_rmaps(:,[lm1,lm2])');

% Compute statistics
[p_value_engaged_vs_moving] = ranksum(engaged_max_dff_at_landmark,moving_max_dff_at_landmark)
[p_value_engaged_vs_still] = ranksum(engaged_max_dff_at_landmark,still_max_dff_at_landmark)
[p_value_moving_vs_still] = ranksum(moving_max_dff_at_landmark,still_max_dff_at_landmark)

engaged_mean = nanmean(engaged_max_dff_at_landmark);
engaged_sem = nanstd(engaged_max_dff_at_landmark)/sqrt(length(engaged_max_dff_at_landmark));

moving_mean = nanmean(moving_max_dff_at_landmark);
moving_sem = nanstd(moving_max_dff_at_landmark)/sqrt(length(moving_max_dff_at_landmark));

still_mean = nanmean(still_max_dff_at_landmark);
still_sem = nanstd(still_max_dff_at_landmark)/sqrt(length(still_max_dff_at_landmark));

bar([engaged_mean,moving_mean,still_mean])
hold on
errorbar([engaged_mean,moving_mean,still_mean],[engaged_sem,moving_sem,still_sem],'LineStyle','none','Color','k')
ylim([0,0.14])
yticks([0,0.14])
ylabel('Mean dF/F at landmark')
xticklabels({'Engaged','Unengaged','Passive'});
xtickangle(25)
set(gca,'FontSize',16);


      
% Set correct size of the figure
set(gcf,'renderer', 'painters', 'Position', [1000,100,500,500])

