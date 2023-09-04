%% FIGURE A:  Plot all cells and highlight the predictive cells
figure(1);
clf;
imagesc(amap_lcs_sorted,[0,3]);
yticks([1,size(amap_lcs_sorted,1)])
ylabel('Cell #')

% Highlight predictive cells
line_width = 1;
line_color = 'r';

% Top
line([0,105],[first_peak_cells(1),first_peak_cells(1)],'Color',line_color,'LineWidth',line_width);
line([0,105],[first_peak_cells(end),first_peak_cells(end)],'Color',line_color,'LineWidth',line_width);
line([1,1],[first_peak_cells(1),first_peak_cells(end)],'Color',line_color,'LineWidth',line_width);
line([105,105],[first_peak_cells(1),first_peak_cells(end)],'Color',line_color,'LineWidth',line_width);

% Bottom
line([0,105],[second_peak_cells(1),second_peak_cells(1)],'Color',line_color,'LineWidth',line_width);
line([0,105],[second_peak_cells(end),second_peak_cells(end)],'Color',line_color,'LineWidth',line_width);
line([1,1],[second_peak_cells(1),second_peak_cells(end)],'Color',line_color,'LineWidth',line_width);
line([105,105],[second_peak_cells(1),second_peak_cells(end)],'Color',line_color,'LineWidth',line_width);

% Landmark onset
line([32,32],[0,2000],'Color','w','LineWidth',line_width)
line([72,72],[0,2000],'Color','w','LineWidth',line_width)
xticks([])

% Set correct size of the figure
set(gca,'FontSize',16)
set(gcf,'renderer', 'painters', 'Position', [2200,100,400,700])