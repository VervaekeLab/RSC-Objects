%% FIGURE D
figure(3);
subplot(1,3,3);
n_predictive_cells = size(predictive_rmaps,3);
n_nonpredictive_cells = size(amap_lcs,1) - n_predictive_cells;


percent = [n_predictive_cells,n_nonpredictive_cells] ./ sum([n_predictive_cells,n_nonpredictive_cells]) * 100;
labels = cellstr( num2str( percent', '%.1f%%' ) );
p = pie([n_predictive_cells,n_nonpredictive_cells], labels);


p(1).FaceColor = [0,0,0];
p(3).FaceColor = [0.8,0.8,0.8];

set(gca,'FontSize',22)



% Set correct size of the figure
set(gcf,'renderer', 'painters', 'Position', [2000,100,800,200])