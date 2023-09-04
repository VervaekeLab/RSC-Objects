%% Parula colorbar
figure(100); clf;
imagesc([0,3],[0,3]);
colormap('parula');
c = colorbar('Ticks',[0,3]);
c.Label.String = "Z-score norm.";
set(gca,'FontSize',16)
set(gcf,'renderer', 'painters', 'Position', [500,100,800,150])

%% Reverse gray colorbar
figure(101); clf;
imagesc([0,3],[0,3]);
colormap(flipud(colormap('gray')));
c = colorbar('Ticks',[0,3]);
c.Label.String = "Z-score norm.";
set(gca,'FontSize',16)
set(gcf,'renderer', 'painters', 'Position', [500,100,800,150])

