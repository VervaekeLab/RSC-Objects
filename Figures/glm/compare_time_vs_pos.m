clearvars -except mData sData data workspace

LLH_increase_pos_mean = [];
LLH_increase_time_mean = [];
pval = [];pval_pos = []; pval_time = [];
for i =1:length(data.sessionIDs)

    pos = load(fullfile(dirct_final,data.sessionIDs{i},'/pos_GLM.mat'));
    time = load(fullfile(dirct_final,data.sessionIDs{i},'/time_GLM.mat'));
    
    for m = 1:length(pos.cell)
        LLH_values_pos = reshape(pos.cell(m).testFit(:,3),5,1);
        LLH_values_time = reshape(time.cell(m).testFit(:,3),5,1);
        LLH_increase_pos_mean(end+1) = mean(LLH_values_pos);
        LLH_increase_time_mean(end+1) = mean(LLH_values_time);
            
          pval_pos(end+1) =signrank(LLH_values_pos,0,'tail','right');
          pval_time(end+1) =signrank(LLH_values_time,0,'tail','right');
          
          if pval_pos(end)<=0.05 & pval_time(end)<=0.05
             pval(end+1) = signrank(LLH_values_pos ,LLH_values_time,'tail','right');
          end
          
    end
       
end


figure()
scatter(LLH_increase_pos_mean,LLH_increase_time_mean,50,'filled','MarkerFaceAlpha',.2)
xlabel('LLH Position')
ylabel('LLH Time')
xlim([-1.5 2])
ylim([-1.5 2])
hold on
plot(-1.5:2,-1.5:2,'k','LineStyle','--','LineWidth',1.5)

ax = gca;
ax.XColor = [0 0 0];
ax.YColor = [0 0 0];
ax.XLabel.Color = [0 0 0];
ax.YLabel.Color = [0 0 0];
set(ax,'FontSize',12,'FontName','Arial')


