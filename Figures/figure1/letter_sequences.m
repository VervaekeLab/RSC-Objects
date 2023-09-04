%% Plot all place cells and landmark cells sorted
amap_pcs = [];
amap_lcs = [];

signal_type = 'deconv_motor_on';
max_z_value = 3;

for s = 1:length(mData)
    
    if s == 1
        amap_pcs = sb.generate.averagedTuningFromRasterMaps(mData(s).rmaps.pcs.(signal_type));
        amap_lcs = [sb.generate.averagedTuningFromRasterMaps(mData(s).rmaps.lcs.(signal_type));sb.generate.averagedTuningFromRasterMaps(mData(s).rmaps.tcs.(signal_type))];
    else
        amap_pcs = cat(1,amap_pcs,sb.generate.averagedTuningFromRasterMaps(mData(s).rmaps.pcs.(signal_type)));
        if ~isempty(mData(s).rmaps.lcs.ind)
            amap_lcs = cat(1,amap_lcs,sb.generate.averagedTuningFromRasterMaps(mData(s).rmaps.lcs.(signal_type)));
        end
        if ~isempty(mData(s).rmaps.tcs.ind)
            amap_lcs = cat(1,amap_lcs,sb.generate.averagedTuningFromRasterMaps(mData(s).rmaps.tcs.(signal_type)));
            
        end
        
        
    end
end

% Sort
amap_pcs = sb.sort.rasterMapByPeak(amap_pcs,'normalize_type','zscore');
amap_lcs = sb.sort.rasterMapByPeak(amap_lcs,'normalize_type','zscore');

% Plot all place and landmark cells sorted

figure()
imagesc([amap_lcs(:,32:end),amap_lcs(:,1:32)],[0,max_z_value]);

hold on
line([0,0],[0,2000],'Color','w');
line([72-32,72-32],[0,2000],'Color','w');

xtix = get(gca, 'XTick');
set(gca, 'XTick',xtix, 'XTickLabel',xtix*1.5)
xticks([1,33,73,105])
xticklabels({'50','100','157','50'})
yticks([1,size(amap_lcs,1)])
set(gca,'FontSize',20);
ylabel('Cell #')
xlabel('Position (cm)')

% Set correct size of the figure
set(gcf,'renderer', 'painters', 'Position', [2500,100,600,500])