
figure()
bar([nanmean(predictive_cells_familiar),...
    nanmean(predictive_cells_introduction),...
    nanmean(precitive_prct_jitter)],'FaceColor','k')
hold on
errorbar(1:3,[nanmean(predictive_cells_familiar),...
    nanmean(predictive_cells_introduction),...
    nanmean(precitive_prct_jitter)],[nanstd(predictive_cells_familiar)/sqrt(length(predictive_cells_familiar)),...
    nanstd(predictive_cells_introduction)/sqrt(length(predictive_cells_introduction)),...
    nanstd(precitive_prct_jitter)/sqrt(length(precitive_prct_jitter))],...
    'LineStyle','none','LineWidth',1.5,'Color',[0.5  0.5 0.5])

hold on
scatter(ones(size(predictive_cells_familiar)),predictive_cells_familiar,50,[0.5  0.5 0.5])
scatter(2*ones(size(predictive_cells_introduction)),predictive_cells_introduction,50,[0.5  0.5 0.5])
scatter(3*ones(size(precitive_prct_jitter)),precitive_prct_jitter,50,[0.5  0.5 0.5])

xticks([1:3])
xticklabels({'Familiar','First object exposure','Unstable object location'})
xtickangle(45)
ylabel('Predictive cells (%)')
set(gca,'FontName','Arial','FontSize',16)
box off
yticks([0 30 60])