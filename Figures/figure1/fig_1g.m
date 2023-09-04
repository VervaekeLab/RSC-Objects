%% This shows the number of cells that are within each class, such as place cells, touch cells etc. 

% Get numbers for all cell categories
n_lcs = sum(data.n_lcs) + sum(data.n_tcs);
n_pcs = sum(data.n_pcs);
n_other = sum(data.n_acs) - (n_lcs+n_pcs);
n_all = n_lcs + n_pcs + n_other;

% Plot
figure(1);
clf;

percent = [n_lcs,n_pcs,n_other] ./ sum([n_lcs,n_pcs,n_other]) * 100;
labels = cellstr( num2str( percent', '%.1f%%' ) );
p = pie([n_lcs,n_pcs,n_other], labels);

% Set color profile
landmark_cells_color = [239,35,60]/255;
place_cells_color = [51,51,51]/255;
other_cells_color = [230,230,230]/255;

p(1).FaceColor = landmark_cells_color;
p(3).FaceColor = place_cells_color;
p(5).FaceColor = other_cells_color;

set(gca,'FontSize',20);

